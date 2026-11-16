@echo off
rem Construye la base de datos de los genes de serie de DNAinvent (~/.dnainvent/catalog.db).
rem Doble clic para los 14 genes. Desde una terminal se le pueden pasar opciones, por ejemplo:
rem     tools\build_catalog.bat --max-variants 200
rem Si se corta, vuelve a ejecutarlo: sigue donde lo dejo.
rem La API key del NCBI (opcional, mas rapido) va en la variable de entorno NCBI_API_KEY,
rem nunca en este archivo.
chcp 65001 >nul
set PYTHONUTF8=1
cd /d "%~dp0.."
python tools\build_catalog.py %*
echo.
pause
