import os
import re
import shutil
import argparse
import textwrap
import subprocess
import importlib.util
import pygments.lexers.python
import pygments.formatters.html
from pathlib import Path
from jinja2 import Environment
from jinja2 import FileSystemLoader
from pdoc.render_helpers import minify_css
import pdoc.render

module_path = Path(__file__).parent.parent
docs_path = module_path / "docs_"
network_viz_path = module_path / "docs" / "network_viz.py"
github_pages_url = "https://isblab.github.io/af_pipeline/"

def get_af_pipeline_version():
    changelog = (module_path / "changelog.md").read_text("utf8")
    # e.g. ## [1.0.0] - 2026/04/06
    version_regex = r"## \[(\d+\.\d+\.\d+)\]"
    for line in changelog.splitlines():
        if line.startswith("## ["):
            match = re.match(version_regex, line)
            if match:
                return match.group(1)

    return "unknown"

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Generate documentation for af_pipeline.")

    parser.add_argument(
        "--output",
        type=str,
        default=str(docs_path),
        help="Output directory for the generated documentation (default: docs_)."
    )

    args = parser.parse_args()
    docs_path = Path(args.output)

    print(get_af_pipeline_version())

    os.makedirs(docs_path, exist_ok=True)
    subprocess.run(
        ["python", str(network_viz_path)],
        check=True
    )

    spec = importlib.util.find_spec("pdoc")
    pdoc_path = Path(spec.origin).parent

    env = Environment(
        loader=FileSystemLoader([module_path, pdoc_path / "templates", pdoc_path / "templates" / "default"]),
        autoescape=True,
    )
    env.filters['minify_css'] = minify_css
    lexer = pygments.lexers.python.PythonLexer()
    formatter = pygments.formatters.html.HtmlFormatter(style="friendly")
    pygments_css = formatter.get_style_defs()

    if docs_path.is_dir():
        shutil.rmtree(docs_path)

    # Render main docs
    pdoc.render.configure(
        edit_url_map={
            "af_pipeline": "https://github.com/isblab/af_pipeline/blob/main/af_pipeline/",
        },
        favicon="./assets/af_pipeline_favicon.svg",
        logo="./assets/af_pipeline_logo.svg",
        logo_link=github_pages_url,
        footer_text=f"af_pipeline v{get_af_pipeline_version()}",
        mermaid=True,
        math=True,
        search=True,
        show_source = True,
        template_directory = module_path / "docs" / "template",
    )

    pdoc.pdoc(
        "af_pipeline",
        output_directory=docs_path,
    )

    # Add sitemap.xml
    with (docs_path / "sitemap.xml").open("w", newline="\n") as f:
        f.write(
            textwrap.dedent(
                """
        <?xml version="1.0" encoding="utf-8"?>
        <urlset xmlns="http://www.sitemaps.org/schemas/sitemap/0.9"
           xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
           xsi:schemaLocation="http://www.sitemaps.org/schemas/sitemap/0.9 http://www.sitemaps.org/schemas/sitemap/0.9/sitemap.xsd">
        """
            ).strip()
        )
        for file in docs_path.glob("**/*.html"):
            if file.name.startswith("_"):
                continue
            filename = str(file.relative_to(docs_path).as_posix()).replace("index.html", "")
            f.write(f"""\n<url><loc>{github_pages_url}{filename}</loc></url>""")
        f.write("""\n</urlset>""")
