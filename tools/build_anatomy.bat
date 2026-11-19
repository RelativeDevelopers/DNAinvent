@echo off
rem Construye «Cuerpo humano» de DNAinvent (~/.dnainvent/anatomy).
rem Doble clic. Necesita los archivos de BodyParts3D en ..\datos\bodyparts3d y, la primera vez,
rem conexion (descarga las tablas de celulas y genes). Opciones, desde una terminal:
rem     tools\build_anatomy.bat --ratio 0.5
chcp 65001 >nul
set PYTHONUTF8=1
cd /d "%~dp0.."
python tools\build_anatomy.py %*
echo.
pause
