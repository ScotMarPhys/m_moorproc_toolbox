# Mooring processing toolbox
Matlab tools for working with instrument data from OSNAP and RAPID-MOC deep water moorings. This toolbox is based on the RAPID-MOC toolbox, the functions contained within this repository take data from raw to gridded and merged time series.  

For documentation,  see Documents/Processing_documents/OSNAP_mooring_processing.pdf

For more info on the data, see [website](https://scotmarphys.github.io/ScotMarPhys.OSNAP-Mooring-Processing.io/)

# How to contribute

### As External
Contribute by opening an issue.

### As User
Before using the toolbox and making any changes to branches it is essential that any user gains a basic understanding of git and github, which can be found here:
https://docs.github.com/en/get-started

# Git version control
This toolbox is managed using git. The master branch should be the last working version used, updated and revised since the last mooring cruise. The master branch is protected and changes can only be merge via pull requests which are then revised by the administrators of the repository. The following strategy is used to manage the repository.

## Updating the master branch with latest cruise branch
1. Before each cruise, users of this mooring processing toolbox should request a new branch from the master branch. Before the users make any changes, they should rename their working branch to reflect the upcoming cruise. E.g. for the cruise JC238 the branch name should be 'JC238'.
2. To keep the master branch up to date, a pull request to merge the cruise branch with the master branch should be raised following good housekeeping practices (below)
   a) Directly after the cruise to push any important on-cruise changes
   b) After on land post-processing when the calibration of the single intruments are finalised and there will be no changes anymore.
3. One of the repository administrators who is not the creator of the branch will review the pull request and check the changes made.
4. If the reviewer finds no issues they will continue to 5. If there are issues or even merge conflicts the reviewer will block the pull request and raise either discussions within the pull request or issues to be addressed by the branch creator. These are open discussion and any contributor are encourage to take part in them to find the best solution. Once an agreement is reached and the creator implemented any nessecary changes and any merge conflicts are resolved steps 2-4 will be repeated.
5. If everything is ok the reviewer can verified the pull request and merge the cruise and master branch.
6. If the cruise branch is finalised and no further changes are expected, it will be locked by one of the administrators making it read only and protecting it from being deleted.

### Working on board
- A copy of the cruise branch should be physically taken on the ship (i.e. osnap/exec/$cruise/) as the internet connection cannot be quarantied on the ship.
- The data (archived from previous cruise) is not kept on git and effort should be made to make sure that the most recent version is copied prior to the cruise/processing. Once the cruise is complete the cruise branch should be pushed to the remote repository (if not already happened during the cruise due to connection problems).

### Housekeeping
It might be nessecary to create feature branches from the cruise branches to solve certain issues (e.g. 'JC238-revise-CM-merging'). When merging the cruise branch into the master branch, any feature branches for the cruise branch should be merged with the cruise branch and deleted so that there is only one branch per cruise. After the final merging (2.b) the cruise branch will be locked by one of the administrators making is read only and no changes can be made.

# Getting started
To set up git and the mooring toolbox please see:

Documents/Processing_documents/OSNAP_mooring_processing.pdf

Please note that the document is work in progress as the mooring toolbox envolves constantly and despite best effort might not always reflect the most current version of the toolbox.

To keep up to date with the most recent updates have also a look on the issues page.
