OptStoic python package
========================
Perform optStoic analysis using Python code that share the same data files with GAMS code.

Note: All the examples are specific for glycolysis pathway generation. 

## Install
- Next, setup a virtual environment in Python 3.
```bash
# Create a project folder
cd project_folder
# Create a virtual environment call optstoic_env
python3 -m venv optstoic_env
# Activate your environment
source optstoic_env/bin/activate
```

- Then, install one of the solvers in the following [Solver Requirement](#solver-requirement) section.

- (Optional) Install the Graphviz package for pathway visualization. See the [Additional Project Dependencies](#additional-project-dependencies) section.

- Next, clone this repository in your `project_folder` and setup. This should install all the Python dependencies.
```
# Create a new project folder
mkdir project_folder
cd project_folder
# Activate your environment
source optstoic_env/bin/activate
# Clone the repo
git clone https://github.com/maranasgroup/optstoic-python.git
cd optstoic-python
python setup.py install
```

- To run nosetests after setup:
```
pip install nose
cd project_folder/optstoic-python
nosetests -s -v
# nosetests -c=nose.cfg
```

## Solver requirement
At least one of the following optimization solvers should be installed. To solve the loopless optStoic formulation, an optimization solver other than GLPK is recommended.

1. GLPK 4.47 installation
   - Linux (Tested on Ubuntu 16.04): 
    ```bash
    wget  http://ftp.gnu.org/gnu/glpk/glpk-4.47.tar.gz
    tar -xvzf glpk-4.47.tar.gz
    cd  ~/glpk-4.47
    ./configure
    make
    make install
    #if the program is successfully installed, you should get an output by typing
    glpsol --version
    ```
    - Mac (Tested on macOS Catalina): 
    ```
    brew install glpk
    # If success
    glpsol --version
    ```

2. GUROBI Optimization provide academic license for free (https://www.gurobi.com/). Install gurobipy following the instruction provided by GUROBI. 

3. [SCIP Optimization Suite](https://scip.zib.de/) >= v4.0.0. See the [documentation of SCIP](https://www.scipopt.org/doc/html/CMAKE.php) for the installation procedure.
    - Linux (Tested on Ubuntu 16.04):
    ```
    sudo apt-get install libgmp-dev libreadline-dev zlib1g-dev libncurses5-dev
    tar xvf scipoptsuite-6.0.0.tgz
    cd scipoptsuite-6.0.0/
    make
    make test
    cd scip-6.0.0/
    sudo make install INSTALLDIR="/usr/local/"
    /usr/local/bin/scip --version
    ```
    - Mac (Tested on macOS Catalina):
    ```
    brew install gmp
    brew install boost
    tar xvf scipoptsuite-7.0.1.tgz
    cd scipoptsuite-7.0.1/
    make
    make test
    cd scip/
    sudo make install INSTALLDIR="/usr/local/"
    /usr/local/bin/scip --version
    ```

4. [CPLEX Optimizer](https://www.ibm.com/analytics/cplex-optimizer)

## Additional project dependencies
1. [PuLP](https://github.com/coin-or/pulp). Run the [test](https://www.coin-or.org/PuLP/main/installing_pulp_at_home.html#testing-your-pulp-installation).

2. Graphviz (Optional, for drawing pathway). The [Graphviz](https://www.graphviz.org/) software is required before installing the graphviz python package. 
    - Linux
    ```bash
    #If you have root access
    sudo apt-get install graphviz

    #If you do not have root access (you can get a different version of Graphviz from their website https://www.graphviz.org/download/)
    cd $HOME
    mkdir -p bin/graphviz
    wget http://www.graphviz.org/pub/graphviz/stable/SOURCES/graphviz-2.38.0.tar.gz
    tar xvf graphviz-2.38.0.tar.gz
    cd graphviz-2.38.0
    ./configure --prefix=$HOME/bin/graphviz
    make && make install
    # Check if the graphviz is working
    cd $HOME/bin/graphviz/bin
    dot -V
    # Add the following line to your .bashrc
    export PATH=$PATH:$HOME/bin/graphviz/bin

    #Install the Python graphviz package
    pip install graphviz
    ```
    - Mac: `brew install graphviz`


3. [Component-Contribution](https://github.com/eladnoor/component-contribution) (*Optional, unless you want to perform MDF analysis)

## Tests
After cloning the repo or setup, please run tests as followed. The runtime depends on the solvers selected by PuLP. Note that the [don't capture stdout](https://nose.readthedocs.io/en/latest/usage.html#cmdoption-s) option must be provided to the nosetests (`nosetests --nocapture` or `nosetests -s`) so that Pulp can read/write from intermediate files.
```
nosetests -s -v
# nosetests --config=nose.cfg
```

## Usage
Read the [tutorial](https://github.com/maranasgroup/optstoic-python/blob/master/optstoicpy/examples/methods.md).

## Jupyter notebook setup
```
cd project_folder
# Activate your environment
source optstoic_env/bin/activate
pip install notebook
pip install ipykernel
python -m ipykernel install --user --name optstoic_env --display-name "Python (optstoic)"
```

## Development
To continue development with the code, please create a virtual environment and use `python setup.py develop` for installation.

## Reference
Please cite [Ng, C.Y., Wang, L., Chowdhury, A. et al. Pareto Optimality Explanation of the Glycolytic Alternatives in Nature. Sci Rep 9, 2633 (2019). https://doi.org/10.1038/s41598-019-38836-9](https://www.nature.com/articles/s41598-019-38836-9).
