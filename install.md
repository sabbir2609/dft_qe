## Installation guide for Quantum Espresso

This guide will help you install Quantum Espresso on your local machine. The guide is written for Ubuntu 22.04, but the steps should be similar for other Linux distributions.

---

### Step 1: **Prepare the Environment**
1. **Update System Packages**:
   ```bash
   sudo apt update && sudo apt upgrade -y
   ```
2. **Install Required Dependencies**:
   - Quantum ESPRESSO requires Fortran and C compilers, libraries for BLAS, LAPACK, and FFTW.
   ```bash
   sudo apt install gfortran gcc make libblas-dev liblapack-dev libfftw3-dev -y
   ```
   Optional: Install MPI for parallel computations:
   ```bash
   sudo apt install mpich libmpich-dev -y
   ```

---

### Step 2: **Download Quantum ESPRESSO**
1. Visit the [official Quantum ESPRESSO website](https://www.quantum-espresso.org/) or the [GitHub repository](https://github.com/QEF/q-e).
2. Download the latest version:
   ```bash
   wget https://github.com/QEF/q-e/archive/refs/tags/qe-7.0.tar.gz
   ```
   7.0 is version number.

3. Extract the tarball:
   ```bash
   tar -xvzf qe-7.0.tar.gz
   cd qe-7.0
   ```

---

### Step 3: **Compile the Source Code**
1. **Configure the Build**:
   Run the configuration script to prepare the makefile:
   ```bash
   ./configure
   ```
   - For parallel execution:
     ```bash
     ./configure MPIF90=mpif90 F90=gfortran
     ```
   If any dependencies are missing, the script will notify you.

2. **Build the Code**:
   Compile the Quantum ESPRESSO suite:
   ```bash
   make all
   ```
   This process might take some time depending on your system.

---

### Step 4: **Set Up Environment Variables**
To use Quantum ESPRESSO executables conveniently, add them to your PATH.

1. Find the installation directory:
   ```bash
   pwd
   ```
   It should be something like `/path/to/qe-X.Y.Z`.

2. Add to PATH:
   ```bash
   echo 'export PATH=/path/to/qe-X.Y.Z/bin:$PATH' >> ~/.bashrc
   ie: echo 'export PATH=/home/sabbir/dft_qe/qe-7.4.1/bin:$PATH' >> ~/.zshrc
   source ~/.bashrc
   ```

---

### Step 5: **Test the Installation**
1. Create a new directory for testing:
   ```bash
   mkdir test_qe
   cd test_qe
   ```
2. ceate a simple input file `input/test.in`: 
   ```bash
   &CONTROL
   calculation = 'scf',
   pseudo_dir = './pseudo',
   outdir = './out',
   /
   &SYSTEM
      ibrav = 2, celldm(1) = 10.2, nat = 2, ntyp = 1,
      ecutwfc = 30.0,
   /
   &ELECTRONS
   /
   ATOMIC_SPECIES
      Si 28.086 Si.pz-vbc.UPF
   ATOMIC_POSITIONS crystal
      Si 0.0 0.0 0.0
      Si 0.25 0.25 0.25
   K_POINTS automatic
      4 4 4 0 0 0
   ```
   - This input file performs a simple SCF calculation for Silicon.
   - pseudo/Si.pz-vbc.UPF is the pseudopotential file for Silicon.
   - out/ is the directory for output files.
   - `Si.pz-vbc.UPF` can be downloaded from the Quantum ESPRESSO website.
   
3. Run the `pw.x` executable:
   ```bash
   pw.x < ./input/test.in > ./output/test.out
   ```
4. Check the output file:
   ```bash
   cat ./output/test.out
   ```
   If the calculation completes without errors, the installation was successful.

---

### Step 6: **Optional - Install Additional Packages**
Quantum ESPRESSO supports various plugins like `EPW`, `Wannier90`, and others. Follow the plugin documentation for installation.

---

### Step 7: **Start Using Quantum ESPRESSO**
- Run `pw.x` for plane-wave calculations:
  ```bash
  pw.x < input_file > output_file
  ```
- Refer to the official documentation for setting up input files and further instructions.

If you encounter any issues during installation, let me know for troubleshooting!