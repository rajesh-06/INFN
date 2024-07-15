#!/bin/bash

# Loop through numbers 00 to 37 (inclusive)
for i in {34..99}; do
  # Create folder name with padding (loop-00, loop-01, ...)
  folder_name="loop-$i"

  # Check if folder already exists
  if [ ! -d "$folder_name" ]; then
    # Create folder if it doesn't exist
    mkdir "$folder_name"
  fi

  # Unzip .tgz file into the folder
  tar -zxvf "loop-$i.tgz" -C "$folder_name"
done

echo "Finished unzipping .tgz files!"
