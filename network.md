
## Network

<a href="../docs/network/af_pipeline_network.html" target="_blank" rel="noopener noreferrer" style="text-decoration: none;">
    <button style="text-align: center; background-color: white; color: #3660a5; border: 2px solid #3660a5; border-radius: 2px; padding: 0.5em; cursor: pointer;" onmouseover="this.style.backgroundColor='#3660a5'; this.style.color='white';" onmouseout="this.style.backgroundColor='white'; this.style.color='#3660a5';">Show in new tab</button>
</a>

<details>

<summary>Click to reveal more information</summary>

<body>
<p>
Double click on the node to go to the corresponding line in the source code.
</p>
<button onclick="openFullscreen()" style="text-align: center; background-color: white; color: #3660a5; border: 2px solid #3660a5; border-radius: 2px; padding: 0.5em; cursor: pointer;" onmouseover="this.style.backgroundColor='#3660a5'; this.style.color='white';" onmouseout="this.style.backgroundColor='white'; this.style.color='#3660a5';">Go Fullscreen</button>

<br>

<iframe id="network_frame" src="../docs/network/af_pipeline_network.html" width="100%" height="800" frameborder="0" allow="fullscreen">
</iframe>

<script>
  function openFullscreen() {
    const elem = document.getElementById("network_frame");
    if (elem.requestFullscreen) {
      elem.requestFullscreen();
    } else if (elem.webkitRequestFullscreen) { /* Safari */
      elem.webkitRequestFullscreen();
    } else if (elem.msRequestFullscreen) { /* IE11 */
      elem.msRequestFullscreen();
    }
  }
</script>
</details>

</body>