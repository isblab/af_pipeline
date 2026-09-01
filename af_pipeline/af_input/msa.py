#!/usr/bin/env python
"""
MIT License

Copyright (c) 2021 Sergey Ovchinnikov

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

import os
import time
import tarfile
import random
import logging
import shutil
from typing import List
import logging
import requests
import argparse
from collections import defaultdict
from tqdm.autonotebook import tqdm
from af_pipeline.constants.af_constants import MMSeqs2API
from af_pipeline.utils.file_utils import read_fasta

logger = logging.getLogger('logger')
TQDM_BAR_FORMAT = '{l_bar}{bar}| {n_fmt}/{total_fmt} [elapsed: {elapsed} remaining: {remaining}]'

class MMseqs2Exception(Exception):
    def __init__(self):
        msg = (
            "MMseqs2 API is giving errors. Please confirm your input is a valid"
            " protein sequence. If error persists, please try again an hour later."
        )
        super().__init__(msg)

def parse_a3m(a3m_file: str) -> dict[str, List[str]]:
    """ Parse .a3m MSA file.

    ## Arguments:

    - **a3m_file (str)**:<br />
        Path to the .a3m file.

    ## Returns:

    - **dict[str, List[str]]**:<br />
        Dictionary where keys are query names and values are lists of aligned
        sequences along with their headers.
    """

    a3m_lines = defaultdict(list)

    with open(a3m_file, "r") as f:
        lines = f.readlines()

    update_query, query = True, None

    for l in lines:

        if len(l.strip()) == 0:
            continue

        if "\x00" in l:
            l = l.replace("\x00", "")
            update_query = True

        if update_query and l.startswith(">"):
            query = l[1:].rstrip()
            update_query = False

        a3m_lines[query].append(l.replace("\t", " "))

    print(a3m_lines.keys())

    return a3m_lines

class MMSeqs2:
    """
    Class to query the MMseqs2 API for MSA generation.
    Adapted from:
        - https://github.com/sokrypton/ColabFold
        - https://github.com/hlasimpk/af3_mmseqs_scripts
    """

    def __init__(
        self,
        sequences,
        targz_file,
        use_env=True,
        use_filter=True,
        use_pairing=False,
        host_url=MMSeqs2API.BASE_URL,
    ):
        self.sequences = sequences
        self.targz_file = os.path.abspath(targz_file)
        self.outdir = os.path.dirname(self.targz_file)
        self.use_env = use_env
        self.use_filter = use_filter
        self.use_pairing = use_pairing
        self.host_url = host_url
        self.setup_mode()

    def setup_mode(self):
        """ Setup the mode for MMseqs2 API based on the provided flags. """

        if self.use_filter:
            self.mode = "env" if self.use_env else "all"
        else:
            self.mode = "env-nofilter" if self.use_env else "nofilter"

        if self.use_pairing:
            self.mode = ""
            self.use_env = False

    def main(self, overwrite=False):
        """ Main function to query MMSeqs2 API, download and extract the results

        ## Arguments:

        - **overwrite (bool, optional):**:<br />
            Whether to overwrite existing results. Defaults to False.
        """

        os.makedirs(self.outdir, exist_ok=True)

        # call mmseqs2 api
        REDO = True

        if os.path.isfile(self.targz_file):
            logger.error(
                f"MMseqs2 results already exist at {self.targz_file}. "
                "Skipping MMseqs2 API call."
            )
            self.extract(targz_file=self.targz_file, overwrite=overwrite)
            return

        # lets do it!
        TIME_ESTIMATE = 150 * len(self.sequences)
        with tqdm(total=TIME_ESTIMATE, bar_format=TQDM_BAR_FORMAT) as pbar:
            while REDO:
                pbar.set_description("SUBMIT")

                # Resubmit job until it goes through
                out = self.submit(self.sequences, self.mode)
                while out["status"] in ["UNKNOWN", "RATELIMIT"]:
                    sleep_time = 5 + random.randint(0, 5)
                    logger.error(f"Sleeping for {sleep_time}s. Reason: {out['status']}")
                    # resubmit
                    time.sleep(sleep_time)
                    out = self.submit(self.sequences, self.mode)

                if out["status"] == "ERROR":
                    raise MMseqs2Exception()

                if out["status"] == "MAINTENANCE":
                    raise MMseqs2Exception()

                # wait for job to finish
                ID, TIME = out["id"], 0
                pbar.set_description(out["status"])
                while out["status"] in ["UNKNOWN", "RUNNING", "PENDING"]:
                    t = 5 + random.randint(0, 5)
                    logger.error(f"Sleeping for {t}s. Reason: {out['status']}")
                    time.sleep(t)
                    out = self.status(ID)
                    pbar.set_description(out["status"])
                    if out["status"] == "RUNNING":
                        TIME += t
                        pbar.update(n=t)

                if out["status"] == "COMPLETE":
                    if TIME < TIME_ESTIMATE:
                        pbar.update(n=(TIME_ESTIMATE - TIME))
                    REDO = False

                if out["status"] == "ERROR":
                    REDO = False
                    raise MMseqs2Exception()

            # Download results
            self.download(ID, self.targz_file)

        self.extract(targz_file=self.targz_file, overwrite=overwrite)

    def submit(
        self,
        sequences: dict[str, str],
        mode: str,
    ) -> dict:
        """ Submit the request to MMSeqs2 API

        ## Arguments:

        - **sequences (dict[str, str])**:<br />
            Dictionary of sequences to be submitted to the MMSeqs2 API.

        - **mode (str)**:<br />
            Mode for the MMSeqs2 API.

        ## Returns:

        - **dict**:<br />
            Dictionary containing the response from the MMSeqs2 API.
        """

        query = ""
        for header, seq in sequences.items():
            query += f">{header}\n{seq}\n"

        _url = MMSeqs2API.get_ticket_url(use_pairing=self.use_pairing)
        res = requests.post(_url, data={"q": query, "mode": mode})

        try:
            out = res.json()

        except ValueError:
            logger.error(f"Server didn't reply with json: {res.text}")
            out = {"status": "ERROR"}

        return out

    def status(self, ID: str) -> dict:
        """ Get status of the submitted request

        ## Arguments:

        - **ID (str)**:<br />
            ID of the submitted request.
            It can be obtained from the output of the `submit` method: out["id"]

        ## Returns:

        - **dict**:<br />
            Dictionary containing the status of the submitted request.
        """

        res = requests.get(MMSeqs2API.get_status_url(ID))

        try:
            out = res.json()
        except ValueError:
            logger.error(f"Server didn't reply with json: {res.text}")
            out = {"status": "ERROR"}

        return out

    def download(self, ID: str, path: str) -> None:
        """ Download the results for the completed request

        ## Arguments:

        - **ID (str)**:<br />
            ID of the completed request.

        - **path (str)**:<br />
            Path to save the downloaded results.
        """

        res = requests.get(MMSeqs2API.get_download_url(ID))

        with open(path, "wb") as out:
            out.write(res.content)

    def extract(self, targz_file: str, overwrite: bool = False) -> None:
        """ Extract the downloaded .tar.gz file to the specified directory.

        ## Arguments:

        - **targz_file (str)**:<br />
            Path to the downloaded .tar.gz file.

        - **overwrite (bool, optional):**:<br />
            Whether to overwrite existing results. Defaults to False.
        """

        tar_name = os.path.basename(targz_file).split(".tar.gz")[0]
        resultdir = os.path.join(os.path.dirname(targz_file), tar_name)
        if os.path.exists(resultdir) and not overwrite:
            logger.error(f"Output directory {resultdir} already exists. Skipping extraction.")
            return

        shutil.rmtree(resultdir, ignore_errors=True)

        with tarfile.open(targz_file, "r:gz") as tar:
            tar.extractall(path=resultdir)


if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--input",
        type=str,
        required=True,
        help="Input sequence file in FASTA format.",
    )

    parser.add_argument(
        "--targz_file",
        type=str,
        required=True,
        help="Output tar.gz file to save the results.",
    )

    parser.add_argument(
        "--use_env",
        action="store_true",
        help="Use environmental sequences.",
    )

    parser.add_argument(
        "--use_filter",
        action="store_true",
        help="Use filtering.",
    )

    parser.add_argument(
        "--use_pairing",
        action="store_true",
        help="Use pairing.",
    )

    args = parser.parse_args()

    # Read input sequences from FASTA file
    sequences = read_fasta(args.input)

    # Create MMSeqs2 object and run the main function
    mmseqs2 = MMSeqs2(
        sequences=sequences,
        targz_file=args.targz_file,
        use_env=bool(args.use_env),
        use_filter=bool(args.use_filter),
        use_pairing=bool(args.use_pairing),
    )
    mmseqs2.main(overwrite=False)