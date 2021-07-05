#!/bin/bash

# generate ssh key
if [ ! -f ~/.ssh/id_ed25519.pub ]; then
    ssh-keygen -t ed25519 -C $emailVar
fi
echo -e "\e[31mssh key =\e[0m $(cat ~/.ssh/id_ed25519.pub)"

# Make sure you have the latest version of the repo
echo
git pull
echo

# # Ask the user for login details
# read -p 'Git repository url: ' upstreamVar
# read -p 'Git Username: ' userVar
# read -p 'Git email: ' emailVar
upstreamVar='git@github.com:willembressers/CarND-Kidnapped-Vehicle-Project.git'
userVar='Willem Bressers'
email='dhr.bressers@gmail.com'

echo
echo Thank you $userVar!, we now have your credentials
echo for upstream $upstreamVar. You must supply your password for each push.
echo

echo setting up git

# specify some global git configuration
git config --global user.name $userVar
git config --global user.email $emailVar
git config --global push.default simple
git config --global core.excludesfile $HOME/.gitignore
git config --global color.ui true
git config --global core.editor vim

# add the remote, specify the branch, pull the latest and greatest
git remote add origin $upstreamVar
git branch --set-upstream-to=origin/main main
git pull origin main
echo

echo Please verify remote:
git remote -v
echo

echo Please verify your credentials:
echo username: `git config user.name`
echo email: `git config user.email`