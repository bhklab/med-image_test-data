#!/bin/bash

# Exit if any command fails
set -e

# Variables
REPO_URL="https://github.com/bhklab/med-image_test-data.git" # Replace with your repository URL
GITHUB_REMOTE="origin" # Remote name, usually 'origin'

# Fetch all tags from the remote
git fetch --tags

# Get the latest tag of the form "v0.##"
LATEST_TAG=$(git tag -l 'v0.*' | sort -V | tail -n 1)

if [ -z "$LATEST_TAG" ]; then
  echo "No tags found in the format 'v0.##'. Starting with v0.01."
  NEW_TAG="v0.01"
else
  # Increment the tag by 0.01
  BASE_VERSION=$(echo "$LATEST_TAG" | cut -c2-) # Remove the 'v' prefix
  NEW_VERSION=$(printf "%.2f" "$(echo "$BASE_VERSION + 0.01" | bc)")
  NEW_TAG="v${NEW_VERSION}"
fi

echo "Latest tag: $LATEST_TAG"
echo "New tag to create: $NEW_TAG"

# Create the new tag
git tag -a "$NEW_TAG" -m "Release $NEW_TAG"

# Push the new tag to the remote repository
git push "$GITHUB_REMOTE" "$NEW_TAG"

echo "Tag $NEW_TAG created and pushed to $GITHUB_REMOTE."