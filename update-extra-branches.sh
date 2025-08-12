#!/bin/bash
# for updating branches that contain a subset of files from another branch.
# useful for projects where you may want to keep the files needed to run something but not various extra features.
## for example, a free version of a paid program, 
## programs with extra files that helped build them but aren't crucial for users, etc.
# while I haven't yet found a way to keep track of dependencies within files, this should suffice.

# WARNING: haven't yet found a way to handle renames where the directory has changed 
## (e.g., copy a directory to a new location and rename it while removing the original), so try to avoid that.

# disable this for now since I don't need it
exit 0

initial_branch=$(git symbolic-ref --short HEAD)
target_branches=("minimal")
if [[ " ${target_branches} " =~ " ${initial_branch} " ]]; then
	echo "Current branch is in the list of branches to update; cancelling script."
	exit 0
fi 
newest_hash=$(git rev-parse HEAD)
files_to_update=$(git diff --name-status --diff-filter=DMRT --find-renames=50% HEAD^ HEAD)
# if no files are updated, exit
if [ -z "$files_to_update" ]; then
	echo "No files to update; cancelling script."
	exit 0
fi
# update each branch
for branch in "${target_branches[@]}"; do
	git checkout $branch
	echo "$files_to_update" | while read -r status file; do
		if [ "$status" == "D" ] || [ "$status" == "M" ] || [ "$status" == "T" ]; then
			# check if file exists
			if [ -f "$file" ]; then
				# if file exists, update it
				git checkout $initial_branch -- $file
			else
				# if file doesn't exist, remove it
				git rm $file
			fi
		else # handle renames
			source_file=$(echo "$file" | cut -f1)
			target_file=$(echo "$file" | cut -f2)
			if [ -f "$source_file" ]; then
				git mv $source_file $target_file
			else
				git rm $target_file
			fi
		fi
	done
	git commit -m "auto update from '$initial_branch' branch"
	git push
done
# return to initial branch
git checkout $initial_branch
