; Inno Setup script: turns dist\MANC-Q into a normal Windows installer (MANC-Q-setup.exe).
; Install Inno Setup (free, jrsoftware.org), open this file, click Compile. Run build_windows.bat first. (GitHub Actions does both automatically.)
#ifndef MyAppVersion
  #define MyAppVersion "1.3.1"
#endif

[Setup]
AppPublisherURL=https://github.com/willwas1/MANC-Q
AppName=MANC-Q
AppVersion={#MyAppVersion}
AppPublisher=University of Manchester
DefaultDirName={autopf}\MANC-Q
DefaultGroupName=MANC-Q
OutputDir=..\dist
OutputBaseFilename=MANC-Q-setup
Compression=lzma2
SolidCompression=yes
PrivilegesRequired=lowest
PrivilegesRequiredOverridesAllowed=dialog
WizardStyle=modern

[Files]
Source: "..\dist\MANC-Q\*"; DestDir: "{app}"; Flags: recursesubdirs createallsubdirs ignoreversion

[Icons]
Name: "{group}\MANC-Q"; Filename: "{app}\MANC-Q.exe"
Name: "{autodesktop}\MANC-Q"; Filename: "{app}\MANC-Q.exe"; Tasks: desktopicon

[Tasks]
Name: "desktopicon"; Description: "Create a desktop shortcut"; Flags: unchecked

[Run]
Filename: "{app}\MANC-Q.exe"; Description: "Start MANC-Q"; Flags: nowait postinstall skipifsilent
