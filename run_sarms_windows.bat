@echo off
echo Installing Julia dependencies...
julia install_dependencies.jl

echo Starting Genie app in production environment...
set GENIE_ENV=prod
julia --project -e "using GenieFramework; Genie.loadapp(); Genie.up();"

echo Opening the app in the default browser...
start http://localhost:8000

echo Genie app is now running.
pause