# EML-electromagnetic-levitation-cell-model
Electromagnetic levitation (EML) cell model for determining the coil and experimental parameters required for stable sample levitation.

A full description of the physics modelled and how it is implemented in the code can be found here: http://hdl.handle.net/2263/61552

## To run an electromagnetic levitation cell simulation using the EML-electromagnetic-levitation-cell-model in Ubuntu

<ul>
  <li>Ensure Python3 and the NumPy, SciPy, and Matplotlib modules for Python3 are installed on the system  </li>
  <pre>
      sudo apt install python3-numpy
      python3 -c "import numpy; print(numpy.__version__)"
  </pre>
  <pre>  
      sudo apt install python3-scipy
      python3 -c "import scipy; print(scipy.__version__)"
  </pre>
  <pre>    
      sudo apt install python3-matplotlib
      python3 -c "import matplotlib; print(matplotlib.__version__)"
  </pre>
  <li>Give the emlc.py file execute permission </li>
  <pre>
    <code class="language-shell">
      chmod +x emlc.py
    </code>
  </pre>
  <li>Run emlc.py with python3 </li>
  <pre>
    <code class="language-shell">
      python3 ./emlc.py
    </code>
  </pre>
</ul>
