# PowerShell script to run run_sarms.jl with Julia and open a browser when "Ready!" is detected

# Determine number of logical processors
$num_cores = (Get-CimInstance Win32_ComputerSystem).NumberOfLogicalProcessors
$env:GENIE_ENV = "prod"

# Check for the presence of installed_libs.so
$libsPath = Join-Path (Get-Location) "installed_libs.so"
if (-Not (Test-Path $libsPath)) {
    # Ask the user if they want to run the installer
    $response = Read-Host "Do you want to install SARMS? This is not essential but will reduce subsequent load times (allow ~10 minutes)."
    if ($response -eq 'Y') {
        Write-Host "Installing SARMS dependencies..."
        $installArgs = "-i", "--threads=$num_cores", ".\install_sarms.jl"
        Start-Process -FilePath "julia" -ArgumentList $installArgs -NoNewWindow -Wait
        Write-Host "Installation complete."
    }
}

Write-Host "Starting SARMS..."

# Prepare the command and arguments
$juliaPath = "julia"  # Make sure Julia is in the PATH or provide the full path to the executable
$arguments = "-i", "--threads=$num_cores", ".\run_sarms.jl"  # Run from the current directory

# Start Julia with the script from the current directory
$processInfo = New-Object System.Diagnostics.ProcessStartInfo
$processInfo.WorkingDirectory = Get-Location  # Set working directory to the current PowerShell path
$processInfo.FileName = $juliaPath
$processInfo.Arguments = $arguments -join " "
$processInfo.RedirectStandardOutput = $true
$processInfo.RedirectStandardError = $true
$processInfo.UseShellExecute = $false
$processInfo.CreateNoWindow = $true

$process = New-Object System.Diagnostics.Process
$process.StartInfo = $processInfo
$process.Start() | Out-Null

# Output handling
$filePath = Join-Path (Get-Location) "sarms_initialiser.log"
$streamWriter = [System.IO.StreamWriter]::new($filePath)
do {
    Start-Sleep -Seconds 1  # Poll every second
    while (!$process.StandardOutput.EndOfStream) {
        $line = $process.StandardOutput.ReadLine()
        $streamWriter.WriteLine($line)
        if ($line -match "Ready!") {
            $ready = $true
            break
        }
    }
} until ($ready -or $process.HasExited)
$streamWriter.Close()

if ($ready) {
    Write-Host "SARMS is now running. Opening browser..."
    Write-Host " ________________________ "
    Write-Host "< SARMS has launched! >"
    Write-Host " ------------------------ "
    Write-Host "        \   ^__^"
    Write-Host "         \  (oo)\_______"
    Write-Host "            (__)       )\/\"
    Write-Host "                ||----w |"
    Write-Host "                ||     ||"

    $browserProcess = Start-Process "http://127.0.0.1:8000" -PassThru
    Write-Host "Monitoring browser process..."

    $filePath = ".\banner.txt"

    # Check if the file exists
    if (Test-Path $filePath) {
        # Read the file line by line
        $asciiArtLines = Get-Content $filePath
        

        
        # Display each line using Write-Host
        foreach ($line in $asciiArtLines) {
            Write-Host $line
            Start-Sleep -Milliseconds 100  # Adjust the delay as needed
        }
    } else {
        Write-Host "File not found: $filePath"
    }

    while (!$browserProcess.HasExited) {
        Start-Sleep -Seconds 1
    }
    Write-Host "Browser has been closed."
    if (!$process.HasExited) {
        $process | Stop-Process -Force
      #  Write-Host "Julia process terminated."

      $juliaProcesses = Get-Process julia -ErrorAction SilentlyContinue

# Check if there are any Julia processes running
if ($juliaProcesses) {
    foreach ($process in $juliaProcesses) {
        Write-Host "Terminating Julia process with PID $($process.Id)..."
        $process | Stop-Process -Force
        Write-Host "Process terminated."
    }
} else {
    Write-Host "No Julia processes found running."
}
    }
} else {
    Write-Host "Failed to detect 'Ready!' in the output. Check logs."
    if (!$process.HasExited) {
        $process | Stop-Process -Force
    }
    exit 1
}

# Wait for user input to close PowerShell window
Write-Host " ________________________ "
Write-Host "< SARMS has closed >"
Write-Host " ------------------------ "
Write-Host "        \   ^__^"
Write-Host "         \  (oo)\_______"
Write-Host "            (__)       )\/\"
Write-Host "                ||----w |"
Write-Host "                ||     ||"


Start-Sleep -Seconds 2

Write-Host " ________________________ "
Write-Host "< Moo >"
Write-Host " ------------------------ "
Write-Host "        \   ^__^"
Write-Host "         \  (oo)\_______"
Write-Host "            (__)       )\/\"
Write-Host "                ||----w |"
Write-Host "                ||     ||"

Write-Host "Press any key to exit ..."


$null = $Host.UI.RawUI.ReadKey("NoEcho,IncludeKeyDown")

exit