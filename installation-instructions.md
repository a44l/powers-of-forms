

## Installation Instructions

Here are some installation instructions, partially written by ChatGPT.  

1. **Install Julia:** Download the latest stable version of Julia from the [official Julia website](https://julialang.org/downloads/). Follow the installation instructions specific to your operating system.

2. **Install Jupyter Notebook:** 
Open a terminal or command prompt. Run the following command to install Jupyter Notebook using Julia's package manager (REPL):
    ```bash
    	julia -e 'using Pkg; Pkg.add("IJulia")'
    ```

3. **Clone the GitHub project:** 
Open a terminal and navigate to the desired directory where you want to clone the project. Run the following command to clone the project:
     ```bash
     git clone https://github.com/a44l/powers-of-forms
     ```

4. **Set up the project environment:**
   - Navigate to the project's directory in the terminal or command prompt.
   - Launch the Julia REPL by running the `julia` command.
   - Type `]` in order to switch to package mode. 
   - Then, type the following commands in the Julia REPL to set up the project environment:

     ```julia
     activate SoSupport
     instantiate
     ```

5. **Start Jupyter Notebook:**
	- In the Julia REPL, exit the package manager mode by pressing the `backspace` key.
	- Type the following command to start Jupyter Notebook:
		```julia-repl
		using IJulia; notebook()
		```
	- A web browser should open automatically, displaying the Jupyter Notebook interface.
    - Navigate to the cloned project directory in the Jupyter Notebook interface.
    - Open `tutorial.ipynb`.
    - Execute the notebook cells to run the code and interact with the project.


### SDP Solver 

Calculations are done using semidefinite programming. Therefore, an interior point solver for SDPs is needed. If no solver is provided, the code will set the `default_solver` to `COSMO`. All experiments and the tutorial notebook explicitly set the solver to `Mosek`. Remember to change the corresponding line, if you do not have access to a `Mosek` licence, but note that free academic licences are available on the [Mosek website](https://www.mosek.com/products/academic-licenses/). The `Mosek` licence must be placed into a folder named `mosek` **in your home directory**. Remember to include the line 
```julia
using MosekTools; default_solver=:Mosek;
``` 
at the beginning of your code files and notebooks, or use the functions `calc_sos_attributes` and `pof_decompose` with the `solver=:Mosek` option. 