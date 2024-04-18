#!/bin/bash

# Script to install dependencies and run a Genie app in production mode

# Get number of cores on macOS
num_cores=$(sysctl -n hw.ncpu)

echo "Installing Julia dependencies..."

echo "Starting SARMS in production environment..."
# Run Julia in interactive mode and attempt to disable output buffering
stdbuf -oL -eL julia -i --threads=$num_cores "run_sarms.jl" &> genie_output.log &
julia_pid=$!

echo "Waiting for SARMS to be ready..."
# Monitor the output log for a specific message indicating readiness
found=0
max_attempts=600
attempts=0

while [ $found -eq 0 ] && [ $attempts -lt $max_attempts ]; do
    if grep -q "Ready!" genie_output.log; then
        found=1
    else
        sleep 1
        ((attempts++))
    fi
done

if [ $found -eq 0 ]; then
    echo "SARMS failed to start within expected time."
    kill $julia_pid
    exit 1
else
    echo "SARMS is now running."
    open http://127.0.0.1:8000
fi

echo " ________________________ "
echo "< SARMS has launched! >"
echo " ------------------------ "
echo "        \\   ^__^"
echo "         \\  (oo)\\_______"
echo "            (__)\\       )\\/\\"
echo "                ||----w |"
echo "                ||     ||"
