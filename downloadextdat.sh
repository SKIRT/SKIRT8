#!/bin/bash
# (use "chmod +rx scriptname" to make script executable)
#
# For use on any Unix system, including Mac OS X and Linux
#
# Execute this script with "git" as default directory to download
# the external resource files provided on the SKIRT web site, and
# place them in the extdat directory next to the git directory.
#

# list of archives to be downloaded and extracted (separate filenames with a space)
FILELIST=( PolarizationProperties.tar.gz FSPSVariableIMFSED.tar.gz )

# select download command: wget (Linux) or curl (Mac OS X)
if which wget >/dev/null
then
    DOWNLOAD="wget --no-check-certificate https://sciences.ugent.be/skirtextdat/SKIRT8/Resources/"
elif which curl >/dev/null
then
    DOWNLOAD="curl --insecure -O https://sciences.ugent.be/skirtextdat/SKIRT8/Resources/"
else
    echo error: no wget or curl available to download files
    exit 1
fi

# download and extract each of the files in the list
filecount=$((${#FILELIST[@]} - 1))
for i in $(eval echo {0..$filecount})
do
    FILENAME=${FILELIST[$i]}
    echo "----------------------------------------"
    echo downloading $FILENAME ...
    mkdir -p ../extdat
    cd ../extdat
    $DOWNLOAD$FILENAME
    tar -xvf $FILENAME
    rm $FILENAME
    cd ../git
done

echo "----------------------------------------"
echo Done.
