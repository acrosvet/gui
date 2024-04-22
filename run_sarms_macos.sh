#!/bin/bash

# Determine number of logical processors
num_cores=$(sysctl -n hw.logicalcpu)
export GENIE_ENV="prod"

# Check for the presence of installed_libs.so
libsPath="./installed_libs.so"
if [ ! -f "$libsPath" ]; then
    read -p "Do you want to install SARMS? (allow ~10 minutes). (Y/n) " response
    if [[ "$response" =~ ^[Yy]$ ]]; then
        echo "Installing SARMS dependencies..."
        if julia -i --threads=$num_cores ./install_sarms.jl; then
            echo "Installation complete. Please restart the script."
            exit 0
        else
            echo "Installation failed."
            exit 1
        fi
    fi
fi


echo "Starting SARMS..."
julia -i --threads=$num_cores ./run_sarms.jl > sarms_initialiser.log 2>&1 &
julia_pid=$!

# Output handling
ready=false
until $ready || ! kill -0 $julia_pid 2>/dev/null; do
    sleep 1
    if grep -q "Ready!" sarms_initialiser.log; then
        ready=true
        break
    fi
done

if $ready; then
    echo "SARMS is now running. Opening browser..."
    echo " ________________________ "
    echo "< SARMS has launched! >"
    echo " ------------------------ "
    echo "        \   ^__^"
    echo "         \  (oo)\_______"
    echo "            (__)       )\/\\"
    echo "                ||----w |"
    echo "                ||     ||"
    
    open "http://127.0.0.1:8000" &
    browser_pid=$!
    echo "Monitoring browser process..."

    filePath="./banner.txt"
    if [ -f "$filePath" ]; then
        while IFS= read -r line; do
            echo "$line"
            sleep 0.1
        done < "$filePath"
    else
        echo "File not found: $filePath"
    fi

    echo "Press any key to close the browser and terminate the Julia process..."
    read -n 1 -s

    if kill -0 $browser_pid 2>/dev/null; then
        echo "Terminating browser process..."
        kill -SIGTERM $browser_pid
        wait $browser_pid 2>/dev/null
        echo "Browser process terminated."
    fi

    if kill -0 $julia_pid 2>/dev/null; then
        echo "Terminating Julia process with PID $julia_pid..."
        kill -SIGTERM $julia_pid
        wait $julia_pid 2>/dev/null
        echo "Julia process terminated."
    fi
else
    echo "Failed to detect 'Ready!' in the output. Check logs."
    if kill -0 $julia_pid 2>/dev/null; then
        echo "Terminating Julia process with PID $julia_pid..."
        kill -SIGTERM $julia_pid
        echo "Julia process terminated."
    fi
    exit 1
fi

echo " ________________________ "
echo "< SARMS has closed >"
echo " ------------------------ "
echo "        \   ^__^"
echo "         \  (oo)\_______"
echo "            (__)       )\/\\"
echo "                ||----w |"
echo "                ||     ||"

echo "Press any key to exit..."
read -n 1 -s
