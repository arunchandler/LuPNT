#!/bin/bash
# 1. Dependencies

echo "Installing dependencies..."
sudo apt-get update
sudo apt-get install -y libboost-all-dev libomp-dev libhdf5-serial-dev
sudo apt install -y qtbase5-dev qtbase5-dev-tools libqt5svg5-dev qttools5-dev-tools

# 2. Download data
echo "Downloading data..."
curl -L https://bit.ly/LuPNT_data -o LuPNT_data.zip
unzip LuPNT_data.zip
rm LuPNT_data.zip

# 3. Set data path based on the shell in use
if [ -n "$ZSH_VERSION" ]; then
    SHELL_RC=~/.zshrc
    SHELL_NAME="zsh"
elif [ -n "$BASH_VERSION" ]; then
    SHELL_RC=~/.bashrc
    SHELL_NAME="bash"
else
    echo "Unsupported shell. Please use bash or zsh."
    exit 1
fi

# Check if LUPNT_DATA_PATH is already set
if [ -n "$LUPNT_DATA_PATH" ]; then
    echo "LUPNT_DATA_PATH is already set to: $LUPNT_DATA_PATH"
    exit 0
fi

# 4. Source the shell configuration file to apply the changes
echo "Applying changes by sourcing $SHELL_RC..."
source $SHELL_RC

# 5. Check data path
echo "LUPNT_DATA_PATH is set to: $LUPNT_DATA_PATH"
