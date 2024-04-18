# PowerShell script to install dependencies and run a Genie app in production mode

# Get number of logical processors
$num_cores = (Get-CimInstance Win32_ComputerSystem).NumberOfLogicalProcessors

Write-Host "Installing Julia dependencies..."

Write-Host "Starting SARMS in production environment..."
# Run Julia in interactive mode and attempt to disable output buffering
$JuliaCmd = "julia -i --threads=$num_cores run_sarms.jl"
Start-Process -FilePath "julia" -ArgumentList "-i", "--threads=$num_cores", "run_sarms.jl" -NoNewWindow -RedirectStandardOutput "genie_output.log" -RedirectStandardError "genie_output.log" -PassThru
$julia_process = Get-Process | Where-Object {$_.ProcessName -eq "julia"}

Write-Host "Waiting for SARMS to be ready..."
# Monitor the output log for a specific message indicating readiness
$found = $false
$max_attempts = 600
$attempts = 0

while (-not $found -and $attempts -lt $max_attempts) {
    Start-Sleep -Seconds 1
    $attempts++
    $logContent = Get-Content "genie_output.log" -Tail 10 -ErrorAction SilentlyContinue
    if ($logContent -match "Ready!") {
        $found = $true
    }
}

if (-not $found) {
    Write-Host "SARMS failed to start within expected time."
    $julia_process | Stop-Process -Force
    exit 1
} else {
    Write-Host "SARMS is now running."
    Start-Process "http://127.0.0.1:8000"
}

Write-Host " ________________________ "
Write-Host "< SARMS has launched! >"
Write-Host " ------------------------ "
Write-Host "        \   ^__^"
Write-Host "         \  (oo)\_______"
Write-Host "            (__)       )\/\"
Write-Host "                ||----w |"
Write-Host "                ||     ||"
