# SEACOAPP

Pre-Installation:
Download the source code (.zip) and unpack it in a folder 
example: Downloads/SEACOAPP-main

Download and unzip the observation data (CMEMS-satellite observation data). 
https://drive.google.com/file/d/13XoX2OxtWl-0VAq4kd9hsRXIuOnOdt-5/view?usp=sharing
Place the unzipped file in the Downloads/SEACOAPP-main.

Installation:
**Windows**
1. Make sure Python 3 is installed
2. Open CMD
3. Write "cd Downloads/SEACOAPP-main" in the cmd prompt. (without quote)
4. Install the requred modules in file reqs.txt
  
  write:
  **pip install -r reqs.txt**
  
  Currently, there is a bug in dash-renderer. 
  Therefore, we recomend to force install older verison of dash-renderer
  
  write again:  
  **pip install --no-dep dash-renderer==1.1.2**
    
3. Run the program
   **python seacoapp.py**
4. The terminal will show the localhost link for the app. Copy the link and paste it on web browser
   example : http://127.0.0.1:XXXX/

**MAC OSX & UBUNTU**
1. Make sure Python 3 is installed
2. Open terminal 
3. Write "**cd Downloads/SEACOAPP-main**" in the terminal. (without quote)
4. Install the requred modules in file reqs.txt
  
  write:
  **pip install -r reqs.txt**
  
  Currently, there is a bug in dash-renderer. 
  Therefore, we recomend to force install older verison of dash-renderer
  
  write again:  
  **pip install --no-dep dash-renderer==1.1.2**
    
3. Run the program
   **python seacoapp.py**
4. The terminal will show the localhost link for the app. Copy the link and paste it on web browser
   example : http://127.0.0.1:XXXX/
