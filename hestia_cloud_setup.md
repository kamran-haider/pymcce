### Steps to configure a jupyter notebook server on Hestia
1. login to hestia

2. I created a shell script to simplify installation of miniconda.
`cp /home/kamran/Dropbox/install_miniconda.sh .`

3. Make sure it is executable.
`chmod +x install_miniconda.sh`

4. Run the script
`./install_miniconda.sh` or `bash install_miniconda.sh`

5. Make sure your .bashrc has a line that points to miniconda python
`vi .bashrc`
add this line
`export PATH="/home/$your_account$/miniconda/bin:$PATH"`
`source .bashrc`

6. Confirm that your python distribution is indeed the one from miniconda
`python`
The header should say Anaconda python 3.6

7. Now we add some useful python packages/libraries. First we add the channels that contains these packages.
```
conda config --add channels conda-forge
conda config --add channels omnia
conda install numpy scipy pandas scikit-learn matplotlib jupyter mdtraj parmed
```

8. Let's also install pymcce
```
conda install pymcce -c kam_haider
```

9. Finally, let's set up a jupyter notebook server that we can use to interact with our files on hestia and create analysis scripts for our data. We have to use different ports for every user
```
jupyter notebook --no-browser --port=8889
```
10. On your local machine, create an ssh tunnel that links localhost on hestia to your localhost.
```
ssh -N -L localhost:8888:localhost:8889 kamran@hestia
```

