@echo off
setlocal
cd /d "%~dp0"
set STREAMLIT_BROWSER_GATHER_USAGE_STATS=false
if not exist ".venv\Scripts\python.exe" (
  echo Please create the environment and install requirements first.
  exit /b 1
)
".venv\Scripts\python.exe" -m streamlit run app.py

