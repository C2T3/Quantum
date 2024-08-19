<h1 align="center">Shor's Algorithm</h1>

### Background
Finding the prime factors of a certain integer $N$ is known to be quite a difficult problem. Indeed, this is such a hard task to do that a widely used encryption system called RSA utilizes this fact to build secure data transmission. By knowing how to factorize numbers efficiently, we could, in theory, break the RSA system which would impact cryptography dramatically. For the moment, there is not a classical algorithm that can factorize any given integer in polynomial time, allowing the RSA encryption system to be robust.

However, with the arrival of the quantum computer, new quantum algorithms exploiting the quantum properties of matter allow for a new way to solve many problems. In fact, some algorithms have been theoretically proven to be more efficient than their classical counterparts. This gives a lot of importance to the development and progress of this new technology.

One of the most promising quantum algorithms is called Shor's algorithm and was developed by MIT professor Peter Shor. It factorizes numbers on a quantum computer with a polynomial complexity, better than any classical algorithm that does the same thing. Since its creation in 1994, the algorithm has gained a lot of attention due to its enormous implications in cryptography as we mentioned earlier.

Still, the implementation of Shor's algorithm seemed at first to be quite challenging because special gates (not trivial ones like the $X$ gate or Hadamard gate for example) were needed for the procedure to work. Many researchers then decided to tackle the issue and published a lot of material to build those special gates. All of this hard work now gives a concrete way to implement Shor's algorithm and to execute it perhaps on real quantum computers.

This repository hosts such an implementation along with extended documentation.

```monospace
├───documentation
|   ├───circuits
|   ├───complements
|   └───shor
|
├───src
|   ├───experiments
|   |   ├───results
|   |   ├───list_semiprimes.npz
|   |   └───semiprimes.py
|   |
|   ├───notebooks
|   |   ├───beauregard.ipynb
|   |   ├───demo.ipynb
|   |   ├───draper.ipynb
|   |   └───pavlidis.ipynb
|   |
|   └───shors_algorithm
|       ├───methods
|       |   ├───utils
|       |   |   ├───draper.py
|       |   |   └───QFT.py
|       |   ├───beauregard.py
|       |   ├───cheat_shor.py
|       |   └───pavlidis.py
|       ├───order_finding.py
|       └───shor.py
|
├───.gitignore
├───LICENSE.txt
├───README.md
└───requirements.txt
```

The folder `documentation` contains many other folders. First, `circuits` has some LaTeX code that was used to build quantum circuits with Quantikz. Finally, `shor` and `complements` both contain LaTeX code that describes the theoretical approach for Shor's algorithm. For the moment, this documentation is only written in French.

The folder `src` contains all the code we have written for Shor's algorithm. There are some notebooks in the folder `notebooks`, some Python files that together implement Shor's algorithm in `shors_algorithm` and material to test the implementation in `experiments`. The header of each Python file has a description of what it does.

### Usage
Create a new virtual environment from `requirements.txt` and feel free to experiment with our implementation after that!

### Credits
This project was done by Mathis Beaudoin (Université de Sherbrooke, Sherbrooke) and Derek Courchesne (C2T3, Trois-Rivières) as a research project for C2T3 during the summer of 2024.

### License
This project is licensed under GPLv3. See `LICENSE.txt` for more details.
