# Developer information

## How to set up the coding environment in Visual Studio Code (Windows)
This section describes how to configure the coding and build toolchain, using the Windows Subsystem for Linux (WSL2), the GCC compiler, CMake and Visual Studio Code:

### 1. Install the required software
- Install the _Windows Subsystem for Linux 2_ (WSL2) as described [here](https://docs.microsoft.com/de-de/windows/wsl/setup/environment) and [here](https://docs.microsoft.com/de-de/windows/wsl/install):
    - Open Windows PowerShell in Admin mode
    - `wsl --install` (installs the default Ubuntu distro)
    - Restart the computer
    - Start the installed WSL2 distribution from the start menu (opens the WSL console)
    - Specify user name and password for the linux admin as prompted
    - `sudo apt update`
    - `sudo apt upgrade`
- Install the _GCC_ compiler and _CMake_ as described [here](https://docs.microsoft.com/de-de/cpp/build/walkthrough-build-debug-wsl2) (the VHoppAS project requires a C++17-compatible compiler and at least CMake version 3.10):
    - `sudo apt install g++ gdb cmake ninja-build rsync zip`
- Download and install _Git for Windows_ (e.g. from [here](https://git-scm.com/download/win))
- Download and install _Visual Studio Code_ on Windows (e.g. from [here](https://code.visualstudio.com/))

### 2. Configure Git (as described [here](https://docs.microsoft.com/de-de/windows/wsl/tutorials/wsl-git))
- Start the _Git Bash_ tool of _Git for Windows_ and enter the following commands to configure Git on Windows (name and email become associated with Git commits; not Github credentials):
    - `git config --global user.name "John Doe"`
    - `git config --global user.email johndoe@example.com`
    - `git config --global credential.helper wincred`
    - `git config --global init.defaultbranch main` (_optional_)
    - `git config --global core.autocrlf true` (_optional_)
- Start the WSL console (Git should be already installed in WSL2) and configure Git in the linux distribution with the same commands, except that the credential helper gets linked to the _Git Credential Helper_ of _Git for Windows_:
    - `git config --global user.name "John Doe"`
    - `git config --global user.email johndoe@example.com`
    - `git config --global credential.helper "/mnt/c/Program Files/Git/mingw64/libexec/git-core/git-credential-wincred.exe"`
    - `git config --global init.defaultbranch main` (_optional_)
    - `git config --global core.autocrlf true` (_optional_)

### 3. Configure Visual Studio Code for coding in WSL2 (as described [here](https://docs.microsoft.com/de-de/windows/wsl/tutorials/wsl-vscode))
- Start the WSL console
- Start _Visual Studio Code_ on Windows
- Install local extensions: `C/C++` and `Remote - WSL`
    - `View` &rarr; `Extensions`
    - Use search bar to find the extensions in the marketplace
    - Install both extensions
- Connect to the WSL linux distribution:
    - `View` &rarr; `Command Palette` &rarr; `Remote-WSL: New WSL Window`
- Install extensions in WSL: `C/C++`, `CMake`, `CMake Tools` and `GitHub Pull Requests and Issues`
    - `View` &rarr; `Extensions`
    - Use search bar to find the extensions in the marketplace
    - Select the extensions and press `Install in WSL`

(_VS Code will re-establish the connection to WSL automatically when it is started while the WSL console is open and when WSL was used in the the previous VS Code session. Otherwise use `Remote-WSL: New WSL Window` in VS Code or enter `code .` in the WSL console._)

### 4. Clone and configure the VHoppAS repository
- In the WSL console create a folder to store repositories, e.g.: `mkdir repos` (_optional_)
- Clone the repository to _Visual Studio Code_ (while connected to WSL): 
    - `View` &rarr; `Command Palette` &rarr; `Git: Clone`
    - Copy & Paste the Github URL of VHoppAS: `https://github.com/phi-hein/VHoppAS.git`
    - Allow redirection to authetication website
    - Enter Github credentials
    - Allow redirection to VS Code in order to transfer the OAuth token
    - Specify the WSL path where to store the repository, e.g.: `/home/`\<_user_\>`/repos/`, and click `OK`
- Configure project-local settings for _Visual Studio Code_:
    - Create a `.vscode` folder in the `VHoppAS` root directory (_this folder name is already included in the `.gitignore` file to exclude it from source control_)
    - Create a `settings.json` file in the `.vscode` folder with at least the following content:
        ```
        {
            "cmake.configureOnOpen": true,
            "C_Cpp.default.cppStandard": "c++17",
            "C_Cpp.default.cStandard": "c17",
            "cmake.buildDirectory": "${workspaceFolder}/build_${buildType}"
        }
        ```
    - Create a `c_cpp_properties.json` file in the `.vscode` folder with at least the following content:
        ```
        {
            "configurations": [
                {
                    "name": "Linux",
                    "includePath": [
                        "${workspaceFolder}/**"
                    ],
                    "defines": [],
                    "compilerPath": "/usr/bin/gcc",
                    "cStandard": "c17",
                    "cppStandard": "c++17",
                    "intelliSenseMode": "linux-gcc-x64",
                    "configurationProvider": "ms-vscode.cmake-tools"
                }
            ],
            "version": 4
        }
        ```
    - Select the configuration: `View` &rarr; `Command Palette` &rarr; `C/C++: Select a Configuration` &rarr; `Linux`
    - Create a `launch.json` file in the `.vscode` folder with at least the following content:
        ```
        {
            "version": "0.2.0",
            "configurations": [
                {
                    "name": "Debug CMake target",
                    "type": "cppdbg",
                    "request": "launch",
                    "program": "${command:cmake.launchTargetPath}",
                    "args": ["-input","DebugInput.txt"],
                    "stopAtEntry": true,
                    "cwd": "${workspaceFolder}/debug",
                    "environment": [
                        {
                            "name": "PATH",
                            "value": "${env:PATH}:${command:cmake.getLaunchTargetDirectory}"
                        }
                    ],
                    "externalConsole": false,
                    "MIMode": "gdb",
                    "setupCommands": [
                        {
                            "description": "Enable pretty-printing for gdb",
                            "text": "-enable-pretty-printing",
                            "ignoreFailures": true
                        }
                    ],
                    "miDebuggerPath": "/usr/bin/gdb"
                }
            ]
        }
        ```
    - Select the _CMake_ compiler kit: `View` &rarr; `Command Palette` &rarr; `CMake: Select a Kit` &rarr; `GCC (General)` (_should refer to the GCC compiler in WSL that supports C++17_)
- Specify VHoppAS input files for debug mode (which define what simulation parameters to use when the VHoppAS is started for example by `Run` &rarr; `Start Debugging`; adapt these files later as needed for debugging):
    - Create a `debug` folder in the `VHoppAS` root directory (_this folder name is already included in the `.gitignore` file to exclude it from source control_).
    - Create a `DebugDOS.txt` file in the `debug` folder, containing a valid DOS specification. For example, copy and rename a DOS file from the `test` folder.
    - Create a `DebugInput.txt` file in the `debug` folder, containing valid simulation parameters. For example, copy and rename an input file from the `test` folder. Set the following file paths in the input file:
        ```
        <MC-Project>
        ...
        DOS-File = "DebugDOS.txt"   
        Output-File = "DebugResult.txt"
        ... 
        </MC-Project>

        ...
        ```
- Choose a build variant: `View` &rarr; `Command Palette` &rarr; `CMake: Select Variant` &rarr; `Debug` or `Release` (_the resulting build folders `build_${buildType}` are already included in the `.gitignore` file to exclude them from source control_)
- Prepare the build files for _CMake_: `View` &rarr; `Command Palette` &rarr; `CMake: Configure` (_this is done automatically for example when the build variant was switched or upon startup of VS Code_)
- Build the VHoppAS program: `View` &rarr; `Command Palette` &rarr; `CMake: Build`
- Run the _CMake_ tests (_optional; may take quite some time to complete_): `View` &rarr; `Command Palette` &rarr; `CMake: Run Tests` (_the resulting output files are already included in the `.gitignore` file to exclude them from source control_)