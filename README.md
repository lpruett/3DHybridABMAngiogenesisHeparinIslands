# 3D Hybrid Agent Based Model of Endothelial Cell Migration in MAP Hydrogel with Heparin Microislands

Lauren Pruett<sup>1</sup>, Alex Taing<sup>1</sup>, Neharika Singh<sup>1</sup>, Shayn Peirce, PhD<sup>1</sup>, Donald Griffin, PhD<sup>1,2,*</sup>

<sup>1</sup>Department of Biomedical Engineering<br>
<sup>2</sup>Department of Chemical Engineering

<sup>*</sup>Corresponding Author<br>
Email address: dg2gf@virginia.edu <br>
Phone number: 434-982-6269<br>

<b>Summary:</b> <br>
A 3D hybrid agent-based model was developed, building off of a previously developed 2D model. The aim of this model
is to simulate endothelial cell migration from a spheroid in a MAP scaffold with heparin microislands in a 3D environment.
The original goal of the 2D model was to optimize the ratio of heparin microislands in a MAP scaffold to promote angiogenesis.
The 3D model aims to do the same, but in a 3D environment.
Heparin microislands in the MAP scaffold were first described here: <br><br>
Pruett, L. J., Jenkins, C. H., Singh, N. S., Catallo, K. J., Griffin, D. R., Heparin Microislands in Microporous
Annealed Particle Scaffolds for Accelerated Diabetic Wound Healing. Adv. Funct. Mater. 2021, 31, 2104337. https://doi.org/10.1002/adfm.202104337

To start using this model, complete installation and setup.

<b>Installation:</b>
1) Install and set up IntelliJ using default settings.
2) Install, set up, and log into GitHub Desktop using default settings.
3) Navigate to http://halloworld.org/index.html in your browser and click "DOWNLOAD". You will be navigated to a GitHub repository.
4) To clone this repository, click on the green "Code" button then select "Open with GitHub Desktop".  You GitHub Desktop will open.  If your browser asks for permissions, allow them.
5) Once GitHub Desktop opens, edit the Local Path to where you would like the library will be installed on your computer. It will be in a file called "HAL-master". _Remember this location for set-up_.
6) Click "Clone". This will install HAL's library on your computer in your specified location.
7) Navigate to this model's repository at https://github.com/lpruett/3DHybridABMAngiogenesisHeparinIslands in your browser.
8) To clone this repository, click on the green "Code" button then select "Open with GitHub Desktop".  You GitHub Desktop will open.  If your browser asks for permissions, allow them.
9) Once GitHub Desktop opens, edit the Local Path to where you would like the library will be installed on your computer. It will be in a file called "3DHybridABMAngiogenesisHeparinIslands". _Remember this location for set-up_.

<b>Set-up:</b>

1) Navigate to the file were the HAL library was cloned to, called "HAL-master". Locate the file where the model was cloned to, called "3DHybridABMAngiogenesisHeparinIslands".
2) Inside "HAL-master", find the folder called "HAL" and the file called "HalColorSchemes.jar". Copy both into the folder called "3DHybridABMAngiogenesisHeparinIslands".
3) Open IntelliJ. If this is your first time opening IntelliJ, press _ctrl + shift + a_ to bring up the search menu. Find "Import Project from Existing Sources" by searching "existing". Then, select the "3DHybridABMAngiogenesisHeparinIslands" folder. Else, if you are currently in another project and not on the IntelliJ "recent projects page", select File then Close Project to be brought back to the project selection window before completing this step.
4) Click accept for default options, <b>stopping when you see "Please select project SDK"</b>. <br> _Note_: If IntelliJ warns about overriding the .idea file, select OK.
5) When you reach "Please select project SDK", click the plus button and download either an SDK with version number 15 or 16.  Once installation finishes, select that new JDK (e.g openjdk-15) in the left window.
6) Click next/create/accept for following options until set-up completes.
7) After loading, you will be greeted by the README. Find the project file tree in the left window.  Select the file "VesselModel3D", then open the file "grid3D".
8) Navigate to "File" in the top left of the screen, followed by "Project Structure", then "libraries".
9) Select the minus button to delete all the current libraries. Click OK on all warnings.
10) Select the plus sign, then select Java. Navigate to the file "lib" located inside the folder "HAL", which itself is inside our "3DHybridABMAngiogenesisHeparinIslands" folder.  After selecting "lib", select OK, then OK on the popup window if one shows.
11) Once again, select the plus sign, then select Java. Navigate to the file called "HalColorSchemes.jar" located directly inside the "3DHybridABMAngiogenesisHeparinIslands" folder. After selecting "HalColorSchemes.jar", select OK.  On the following popup box, select "Jar Directory", then OK. <br>
    _Note_: this may seem like you just added the libraries you removed, but these steps eliminate an error with programs that use the OpenGL3DWindow class.
12) Select OK on the "Project Structure" window.
13) Inside the "3DHybridABMAngiogenesisHeparinIslands" folder, navigate to HAL/GridsAndAgents/SphericalAgent3D.java.  Additionally, locate the Required_additions.txt file located inside VesselModel3D. Inside SphericalAgent3D.java, add all text below the dashed line (line 4 and below) from Required_additions.txt below the fields of the class (to line 30, which should be empty).  You may save and close this window.
14) Click "Add Configuration" located in the top right corner between the hammer button and the play button.
15) In the opened window, click the plus sign in the upper left corner to open a dropdown menu. Select application.
16) Name the application (suggested: VesselModel3D).
17) Under the "Build and Run" header, select the main class called "grid3D" by clicking on the file icon at the end of the first input box and searching for "grid3D".
18) Select OK. Now the program can be run by selecting the play button in the top right corner. <br>
    _Note_: If an error occurs, carefully complete steps 10 and 11 once more and try running the program again.

<br>

<b>Acknowledgements:</b> <br><br>
HAL: http://halloworld.org/index.html <br>
Bravo RR, Baratchart E, West J, Schenck RO, Miller AK, Gallaher J, et al. (2020) Hybrid Automata Library: A flexible platform for hybrid modeling with real-time visualization. PLoS Comput Biol 16(3): e1007635. https://doi.org/10.1371/journal.pcbi.1007635
<br>
<br>
<br>
This project was supported by the UVA Double Hoo grant. <br>
L. Pruett is funded by NIH F31HL154731<br>
This work was partially supported by NIH R01 10297936<br>
