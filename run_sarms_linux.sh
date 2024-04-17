#!/bin/bash

# Script to install dependencies and run a Genie app in production mode

num_cores=$(nproc)

echo "Installing Julia dependencies..."

echo "Starting SARMS in production environment..."
# Run Julia in interactive mode and attempt to disable output buffering
stdbuf -oL -eL julia -i --threads=$num_cores "run_sarms.jl" &> genie_output.log &
julia_pid=$!

echo "Waiting for SARMS to be ready..."
# Monitor the output log for a specific message indicating readiness
while ! grep -q "Ready!" genie_output.log; do
    sleep 1
done

echo "SARMS is now running."
xdg-open http://127.0.0.1:8000 &

echo " ________________________ "
echo "< SARMS has launched! >"
echo " ------------------------ "
echo "        \\   ^__^"
echo "         \\  (oo)\\_______"
echo "            (__)\\       )\\/\\"
echo "                ||----w |"
echo "                ||     ||"
