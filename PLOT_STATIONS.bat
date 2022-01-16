cd C:\EQP11-20\Programs
call C:\Users\DELL\anaconda3\Scripts\activate.bat pygmt
python PLOT_STATIONS_PYTHON.py
call conda deactivate
python pdf_merge.py
cd C:\Program Files (x86)\Google\Chrome\Application
call chrome.exe C:\EQP11-20\Programs\Result.pdf
cd C:\EQP11-20\Programs

