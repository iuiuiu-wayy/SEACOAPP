# SEACOAPP

## Southeast Asia Climate-Ocean Web App 

### Developed by:
 * Imam Wahyu Amanullah, Universiti Teknikal Melaka Malaysia, Malaysia
 * Fadhlil Rizki Muhammad, The University of Melbourne, Australia
 * Akhmad Faqih, IPB University, Indonesia

### **BAHASA INDONESIA**

Pra-Instalasi:
Unduh kode sumber (.zip) dan ekstrak dalam sebuah folder
contoh: Downloads / SEACOAPP-main

Unduh dan unzip data observasi (CMEMS-data observasi satelit).
https://drive.google.com/file/d/13XoX2OxtWl-0VAq4kd9hsRXIuOnOdt-5/view?usp=sharing
Tempatkan file unzip di Downloads / SEACOAPP-main.

Instalasi:

**Windows**
1. Pastikan Python 3 diinstal
2. Buka CMD
3. ketik "cd Downloads / SEACOAPP-main" pada prompt cmd. (tanpa tanda kutip)
4. Instal modul yang diperlukan di file reqs.txt
  
  ketik:
  **pip install -r reqs.txt**
  
  Saat ini, ada bug di dash-renderer.
  Oleh karena itu, kami menyarankan untuk memasang paksa versi lama dari dash-renderer
  
  tulis lagi:
  **pip install --no-dep dash-renderer==1.1.2**
    
3. Jalankan programnya
   **python seacoapp.py**
4. Terminal akan menampilkan link localhost untuk aplikasi tersebut. Salin tautan dan tempel di browser web
   
   contoh: http://127.0.0.1:XXXX/

**MAC OSX & UBUNTU**

1. Pastikan Python 3 diinstal
2. Buka terminal
3. Tulis "**cd Downloads / SEACOAPP-main**" di terminal. (tanpa kutipan)
4. Instal modul yang diperlukan di file reqs.txt
  
  ketik:
  **pip install -r reqs.txt**
  
  Saat ini, ada bug di dash-renderer.
  Oleh karena itu, kami menyarankan untuk memasang paksa versi lama dari penyaji dasbor
  
  tulis lagi:
  **pip install --no-dep dash-renderer == 1.1.2**
    
3. Jalankan programnya
   **python seacoapp.py**
4. Terminal akan menampilkan tautan localhost untuk aplikasi tersebut. Salin tautan dan tempel di browser web
   
   contoh: http://127.0.0.1:XXXX/

### **ENGLISH**

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
