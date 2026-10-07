@echo off
REM ============================================================================================
REM  MANC-Q launcher for Windows. Double-click this file to start MANC-Q.
REM  From a Command Prompt, "MANC-Q.bat run ..." etc. runs the command line instead (see docs\COMMAND_LINE.md).
REM
REM  The first time, it sets everything up inside this folder (about 5 minutes, needs internet):
REM    1. downloads "uv", a small installer tool, into MANC-Q\.runtime
REM    2. uv downloads a private copy of Python and the packages MANC-Q needs, also into .runtime
REM  Nothing is installed on the rest of the computer and no administrator rights are needed.
REM  After that it opens in a few seconds. To uninstall, delete the MANC-Q folder.
REM ============================================================================================
setlocal
title MANC-Q (keep this window open while MANC-Q is running)
set "CALLDIR=%CD%"
cd /d "%~dp0"

set "RT=%~dp0.runtime"
set "VENV=%RT%\venv"
set "UV=%RT%\uv\uv.exe"
set "UV_CACHE_DIR=%RT%\cache"
set "UV_PYTHON_INSTALL_DIR=%RT%\python"
set "UV_PYTHON_PREFERENCE=only-managed"
set "UV_NATIVE_TLS=1"

if not exist "pyproject.toml" (
    echo This file must stay inside the MANC-Q folder, next to pyproject.toml.
    echo If you opened it from inside the ZIP file, extract the ZIP first ^(right-click, Extract All^).
    goto :fail
)

REM Already set up, and the MANC-Q files have not changed since? Then just start.
if exist "%VENV%\Scripts\python.exe" if exist "%RT%\installed_stamp.txt" (
    copy /b "pyproject.toml"+"constraints.txt" "%RT%\now.txt" >nul && fc /b "%RT%\now.txt" "%RT%\installed_stamp.txt" >nul 2>&1 && goto :start
)

echo.
echo  Setting up MANC-Q in this folder. This happens once and takes about 5 minutes.
echo  ------------------------------------------------------------------------------
if not exist "%RT%" mkdir "%RT%"

if not exist "%UV%" (
    echo  [1/3] Downloading the installer tool ^(uv^)...
    if not exist "%RT%\uv" mkdir "%RT%\uv"
    curl.exe -L --fail --silent --show-error -o "%RT%\uv.zip" https://github.com/astral-sh/uv/releases/latest/download/uv-x86_64-pc-windows-msvc.zip || goto :nonet
    tar -xf "%RT%\uv.zip" -C "%RT%\uv" || goto :fail
    del "%RT%\uv.zip"
)

echo  [2/3] Getting Python ^(a private copy, only used by MANC-Q^)...
if exist "%VENV%" rmdir /s /q "%VENV%"
"%UV%" venv --python 3.12 "%VENV%" || goto :nonet

echo  [3/3] Installing MANC-Q and the packages it uses...
"%UV%" pip install --python "%VENV%\Scripts\python.exe" -e . -c constraints.txt || goto :nonet
copy /b "pyproject.toml"+"constraints.txt" "%RT%\installed_stamp.txt" >nul

REM A shortcut on the desktop, so next time you do not need to find this folder.
powershell -NoProfile -Command "$s=(New-Object -ComObject WScript.Shell).CreateShortcut([Environment]::GetFolderPath('Desktop')+'\MANC-Q.lnk'); $s.TargetPath='%~f0'; $s.WorkingDirectory='%~dp0'; $s.Save()" >nul 2>&1 && echo  A MANC-Q shortcut has been put on your desktop.
echo  Set-up finished.
echo.

:start
if not "%~1"=="" goto :command
echo  Starting MANC-Q. The window can take up to a minute to appear the first time.
echo  Leave this black window open; closing it closes MANC-Q.
"%VENV%\Scripts\python.exe" -m mancq gui
if errorlevel 1 goto :fail
exit /b 0

:command
REM  Anything typed after MANC-Q.bat is passed to the command line, e.g.
REM    MANC-Q.bat run "C:\NMR\Experiment 1" "C:\NMR\Results" --ref-mm 0.5
cd /d "%CALLDIR%"
"%VENV%\Scripts\python.exe" -m mancq %*
exit /b %errorlevel%

:nonet
echo.
echo  The download failed. Check that this computer is connected to the internet. On a university or
echo  company network the firewall may block downloads; try another network once, or use the installer
echo  ^(MANC-Q-setup.exe^) from the Releases page: https://github.com/willwas1/MANC-Q/releases
echo  Run this file again to retry; it carries on where it stopped.
goto :fail

:fail
echo.
echo  MANC-Q stopped because of the problem above. If you ask for help, send a screenshot of this window.
pause
exit /b 1
