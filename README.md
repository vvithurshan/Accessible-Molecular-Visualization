<div style="display: flex;">
  <img src="https://github.com/vvithurshan/Pymol/blob/main/images/pymol.jpg" width="50" height="50">
  <img src="https://github.com/vvithurshan/Pymol/blob/main/images/windows.png" width="50" height="50">
  <img src="https://github.com/vvithurshan/Pymol/blob/main/images/mac.jpeg" width="50" height="50">
  <img src="https://github.com/vvithurshan/Pymol/blob/main/images/ubuntu-logo-clipart-5.png" width="50" height="50">
</div>

 <!--
<div style="display: flex;">
  <div style="margin-right: 10px;">
    <img src="https://github.com/vvithurshan/Pymol/blob/main/images/pymol.jpg" width="50" height="50">
  </div>
  <div style="margin-right: 10px;">
    <img src="https://github.com/vvithurshan/Pymol/blob/main/images/windows.png" width="50" height="50">
  </div>
  <div style="margin-right: 10px;">
    <img src="https://github.com/vvithurshan/Pymol/blob/main/images/mac.jpeg" width="50" height="50">
  </div>
  <div>
    <img src="https://github.com/vvithurshan/Pymol/blob/main/images/ubuntu-logo-clipart-5.png" width="50" height="50">
  </div>
</div>
src="https://github.com/vvithurshan/Pymol/blob/main/images/ubuntu-logo-clipart-5.png" width="50" height="50">
</div>
-->




# How to use 
## 1. Cloning the Repository
Open a terminal or command prompt.

Navigate to the directory where you want to clone the repository.

Use the git clone command followed by the repository URL:
```bash
  git clone https://github.com/vvithurshan/Pymol.git
```

Press Enter to execute the command.

Wait for Git to clone the repository to your local machine.

Once the cloning process is complete, you will have a local copy of the repository in your specified directory.

## 2.  Install Dependencies

To install the required packages, execute the following command where you cloned the repository:
    
**For Windows OS**
```bash
pip install -r requirments.txt
```
**For Linux and Mac OS**
```bash
pip3 install -r requirments.txt
```

## 3. Using PyMOL
 #### Follow these steps to set up the environment:

    1. Open PyMOL.
    2. Go to the "File" menu.
    3. Select "Run Script".
    4. Choose "pymol.pymol" from the file dialog.
    5. Repeat steps 2-4, but this time select "server.py" from the file dialog.
    6. This will load the necessary scripts for the application.

## 4. In order to use the large language model called gemini-pro, API key has to be created in the following link.
	```http
 	https://aistudio.google.com/app/prompts/new_chat
 	```
  ### After getting the API key, the variable called GOOGLE_API_KEY (line number 395) has to be updated in the file called application.py
## 5. Now, open a terminal and navigate to the directory where the repository is cloned. Then, start the Flask application by executing the following command:
**For Windows OS**
```bash
python application.py
```
**For Linux and Mac OS**
```bash
python3 application.py
```
**If the Flask application starts successfully, you should see a similar message in the terminal:**

```bash
 * Serving Flask app 'application'
 * Debug mode: on
WARNING: This is a development server. Do not use it in a production deployment. Use a production WSGI server instead.
 * Running on all addresses (0.0.0.0)
 * Running on http://127.0.0.1:5000
 * Running on http://10.132.10.89:5000
Press CTRL+C to quit
 * Restarting with stat
 * Debugger is active!
 * Debugger PIN: 144-667-037
```
## 6. After launching the Flask application, open a web browser and navigate to the following URL in the address bar: 

```
http://localhost:5000/
```
### The following menu appears in the web browser if successful
<div>
<img src="https://github.com/vvithurshan/Pymol/blob/main/images/front_menu.png" width="600" height="400">
</div>

## 7. We can extend the package's functionality by building upon it. Here are some ways to achieve this:
- The available functionality available in pymol can be known from the following web page.
  ```
  https://pymol.org/pymol-command-ref.html
  ```

- If you want to add any functionality, you can simply write a function in the file called pymol.py as follows.
  ```python
  ## rotate in x axis
  def rotate_x(deg):
	  cmd.rotate('x', deg)

  ## Enable Depth Cue
  def Depth_Cue():
	  cmd.set("depth_cue", 1)
	  cmd.set("fog_start", 0.03)
  ```

- After writing the functions in pymol.py, the name of the functions have to be mentioned in the file called server.py as shown below.
  ```python
  httpd.expose("rotate_x", rotate_x)
  httpd.expose("Depth_Cue", Depth_Cue)
  ```
- The application.py file serves as the core of our web application's interaction with the Pymol server.  Here, we'll utilize a conditional statement (elif) within the output() function to handle user input effectively. 

```python 
    #Zoomin
    elif input_data == 'zoomin' or input_data == 'zoom in':
        if pdb_code:
            respose = requests.post('http://localhost:8080/apply/zoomin')
            return render_template('index.html', final_result = "Zooming in")

```
- Capturing User Input: The input_data variable plays a crucial role. It captures the user's input submitted through the web application. This data will determine the specific action taken by the code.

- Conditional Logic (elif): The elif statement provides a way to execute specific code blocks based on the value of input_data.  Here's a breakdown of the logic:

- The main output() function likely has an initial if statement handling a specific case.The elif  statement provides an alternative scenario based on the value of input_data. Inside the elif block, you'll write the URL for the HTTP request that needs to be sent to the Pymol server. This URL will likely be constructed dynamically based on the user's input.

- Rendering the Output: Once the elif block executes the HTTP request and retrieves a response from the Pymol server, you'll likely want to display the results for the user. Here's how the render_template function comes into play:

- The first argument, 'index.html', specifies the HTML template file that will be used to render the final output. This template likely contains placeholders for dynamic content. The second argument, final_result, holds the data retrieved from the Pymol server. This data will be passed to the template and displayed in the designated text box within the index.html file.

**This approach allows the web application to interact with the Pymol server dynamically based on user input.  The results are then displayed visually within the webpage and potentially provided as voice output for an enhanced user experience.**

## Available Commands
1. load - for loading pdb files.
```
## for loading a pdb file called 4cha
load 4cha
```

2. zoomin

3. zoomout

4. rotatex
5. rotatey
6. rotatez
7. reset
8. marker
9. delete
10. reset
11. rg
12. secondary
13. contact map
14. llm
