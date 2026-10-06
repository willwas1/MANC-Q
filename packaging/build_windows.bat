@echo off
REM Builds the stand-alone MANC-Q application into dist\MANC-Q\MANC-Q.exe on your own PC.
REM You only need this if you want to build it yourself: GitHub builds it automatically when you publish a
REM release (see .github\workflows\build-windows.yml). Needs Python 3.9+ on the PATH.
REM Run from the MANC-Q folder:   packaging\build_windows.bat

python -m venv build_env || goto :error
call build_env\Scripts\activate.bat || goto :error
python -m pip install --upgrade pip || goto :error
python -m pip install . pyinstaller -c constraints.txt || goto :error
python packaging\build_app.py || goto :error

echo.
echo Self-test: fitting the synthetic demo spectra with the built application (a few minutes)...
start /wait "" dist\MANC-Q\MANC-Q.exe --selftest
type "%TEMP%\MANC-Q selftest\selftest.txt"

echo.
echo Done. The application is in dist\MANC-Q\  (start dist\MANC-Q\MANC-Q.exe)
echo To make an installer, open packaging\installer.iss in Inno Setup and click Compile.
call deactivate
exit /b 0

:error
echo.
echo The build failed at the step above.
exit /b 1
