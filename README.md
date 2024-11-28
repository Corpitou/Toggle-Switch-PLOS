# Toggle Switch Model 

## Description
This project implements the Gillespie Stochastic Simulation Algorithm to model the dynamics of a two-gene regulatory network. The focus is on the robustness of differentiation processes and the influence of noise within a genetic regulatory system. The model encompasses multiple reaction pathways, represented by rate equations that depend on parameters and variables describing gene activation and degradation.

A XPP-AUTO file is also implemented to compute Temporal Series, Phase Spaces and Bifurcation Diagrams (see [XPP-AUTO Website](https://sites.pitt.edu/~phase/bard/bardware/xpp/xpp.html) for more informations.)

---

## Table of Content
- [Description](#description)
- [Quickstart for Advanced Users](#quickstart-for-advanced-users)
- [Full Installation Guide](#full-installation-guide)
  - [Prerequisites](#prerequisites)
  - [Platform-Specific Instructions](#platform-specific-instructions)
    - [macOS](#macos)
    - [Linux](#linux)
    - [Windows (via WSL)](#windows-via-wsl)
    - [Notes](#notes)
- [Usage](#usage)
  - [Compilation and running](#compilation-and-running)
  - [Customizing Programs](#customizing-programs)
- [Troubleshooting](#troubleshooting)
  - [Common Problems](#common-problems)
    - [macos Users](#macos-users)
    - [Linux Users](#linux-users)
    - [Windows Users](#windows-users)
    - [Useful Ressources](#useful-ressources)
- [GSL Installation](#gsl-installation)
  - [macOS](#macos)
    - [(Recommended) Website Installation](#recommended-website-installation)
    - [(Optional) Homebrew Installation](#optional-homebrew-installation)
  - [Linux](#linux)
- [General Note](#general-note)
- [Contact](#contact)
- [License](#license)

---

## Quickstart for Advanced Users

1. **Ensure prerequisites are installed:**
   - **Compiler and tools:** Ensure that `gcc`, `make`, and basic development tools are installed on your system.
   - **GNU Scientific Library (GSL):** Although not strictly required, having GSL installed on your system (e.g., `gsl-config --version`) can help with debugging. Precompiled GSL static libraries are already included in the `LIB` directory of the project.

2. **Clone the repository and compile:**
   ```bash
   git clone https://github.com/Corpitou/Toggle-Switch-PLOS.git
   cd Toggle-Switch-PLOS
   make

3. **Compilation and execution:**

   - The make command will compile the project and automatically run the program. If the compilation is successful, you will see the output: `Hello, it's working!`.

3. **Modify the code:**

   - Open SRC/main.c, make changes to the main() function or other files, and recompile using make.

---

## Full Installation Guide

### Prerequisites
- **Operating System:** macOS (M-chip recommended), Linux, or Windows (via WSL).
- **Required Tools:**
  - GNU Scientific Library (GSL) (Although not strictly required, having GSL installed on your system (e.g., `gsl-config --version`) can help with debugging. Precompiled GSL static libraries are already included in the `LIB` directory of the project.)
  - C compiler (gcc or clang)
  - Make utility
- **Optional:**
  - [Homebrew](https://brew.sh) (not recommended for GCC or GSL due to potential conflicts; prefer downloading GCC from the XCode Command Line Tools and downloading GSL directly from the official site).
  
### Platform-Specific Instructions

#### macOS
1. **(Recommended) Install Xcode Command Line Tools:**

The easiest way to get gcc or clang on macOS is by installing the Command Line Tools for Xcode:

   ```bash
   xcode-select --install
   ```
This will install 'clang', 'make', and other essential tools.

2. **(Optional) Install Homebrew:**

Homebrew is a package manager for macOS and can be used to install a variety of software, including gcc and other development tools. To install Homebrew:

   ```bash
   /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
   ```
   Use Homebrew to install dependencies:
   ```bash
   brew install gcc make gsl
   ```
3. **Verify installations:**

   ```bash
   gcc --version
   make --version
   gsl-config --version
   ```

#### Linux
On Linux, the required tools and libraries can be installed using your system's package manager. Below are the instructions for systems based on Debian/Ubuntu.

1. **Update your package manager:**
   ```bash
   sudo apt update
   ```
2. **Install required tools:**
   ```bash
   sudo apt install -y gcc make build-essential libgsl-dev
   ```
3. **Verify installations:**
   ```bash
   gcc --version
   make --version
   gsl-config --version
   ```

#### Windows (via WSL)
1. **Install WSL:**
   ```bash
   wsl --install
   ```
   Restart your system if prompted.
2. **Open a Linux terminal via WSL and install dependencies:**
   ```bash
   sudo apt update
   sudo apt install libgsl-dev build-essential
   ```
3. **Clone the repository and follow Linux instructions.**

#### Notes
macOS vs Linux Compatibility: The project has been tested and is compatible with both macOS (using gcc or clang) and Linux (gcc).
Static Libraries: The project uses pre-compiled static libraries for GSL, so you don't need to install dynamic GSL libraries. However, installing GSL on your system may help with debugging or future modifications.

---

## Usage

### Compilation and Running

1. **Open a terminal.**
2. **To compile and run this project, navigate to the project directory and use the provided Makefile:**
```bash
cd /path/to/Toggle-Switch-PLOS
make
```
This will execute the script in the Makefile, detect the GSL installation path, and compile the project. If everything works correctly, the compiled program will automatically run and display "Hello, it's working!".

### Customizing Programs

1. **Open the `main.c` file (in `SRC/`).**
2. **Modify the `main()` function to call the desired function (e.g., `START_...`).**
3. **Recompile the code using `make`.**

---

## Troubleshooting

### Common Problems

#### macOS Users

Normally, the project should run properly on macOS (at least with a M-chip). 

#### Linux Users

Some issues may arise on Linux, particularly with **basic includes**. Normally, the basic libraries that are called in `main_lib.h` are selected depending on the OS (macOS or Linux). 

If you encounter compilation or runtime errors, check:
- That the **GNU Scientific Library (GSL)** is properly installed.
- That all standard includes (e.g., `stdlib.h`, `stdio.h`) are supported on your system.
- Ensure that your compiler supports the required flags and libraries.

#### Windows Users
Use WSL (Windows Subsystem for Linux) to run the project. Follow the [WSL instructions](#windows-via-wsl) above.

#### Useful Ressources
- Official GSL Documentation: [https://www.gnu.org/software/gsl/](https://www.gnu.org/software/gsl/)
- Stack Overflow for specific issues: [https://stackoverflow.com/](https://stackoverflow.com/)

---

## GNU Scientific Library (GSL) Installation

### macOS

#### (Recommended) Website installation 

1. **Directly install the GNU Scientific Libraries from the website: [https://www.gnu.org/software/gsl/](https://www.gnu.org/software/gsl/)**

2. **Follow the `INSTALL` file instructions.**

#### (Optional) Homebrew installation 

1. **Install or reinstall GSL using Homebrew:**
   ```bash
	brew install gsl
   ```

2. **Verify the installation:**
   ```bash
	gsl-config --version
   ```

As explain before, it is not recommended for GSL due to potential conflicts; prefer downloading GSL directly from the official site.

### Linux
1. **Install GSL using your package manager:**

   ```bash
	sudo apt update
	sudo apt install libgsl-dev
   ```

---

## General Note

If you encounter problems, please:

1. Check the documentation and links above.
2. Note that **this project will not receive future support or updates.**

---

## Contact
For questions or issues, contact:
**Email:** corentin.robert@ulb.be

---

## License
This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

---

