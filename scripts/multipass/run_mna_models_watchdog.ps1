
Start-Job -ScriptBlock {
    while ($true) {
        Get-Process | Where-Object { $_.MainWindowTitle -match "Ошибка|Error|Crash" } | ForEach-Object { Stop-Process -Id $_.Id -Force }
        Start-Sleep -Seconds 1
    }
} | Out-Null

# Define paths
$exePath = ".\OLMPASSdoMNABase.exe"
$inputFolder = ".\uncharged"
$outputBaseFolder = "MSAR"  # Base folder for MSAR directories (assumes they already exist)

# Clear or create the consolidated output file
$consolidatedOutputFile = ".\ConsolidatedOutput.txt"
New-Item -Path $consolidatedOutputFile -ItemType File -Force | Out-Null

# Get all .sdf files in the input folder
$sdfFiles = Get-ChildItem -Path $inputFolder -Filter "*.sdf"

foreach ($file in $sdfFiles) {
    # Extract protein name and number from the file name
    $fileName = $file.Name
    if ($fileName -match "^(?<protein>\w+)_(?<number>\d+)\.sdf$") {
        $protein = $matches['protein']
        $number = $matches['number']

        # Loop through numbers 3 to 30 for the BaseCreate command
        for ($i = 3; $i -le 30; $i++) {
            # Define the MSAR directory for the current number
            $msarDirectory = ".\MSAR_$i"

            # Check if the MSAR directory exists, then proceed
            if (Test-Path -Path $msarDirectory) {
                # Create a DoListTask file for the current iteration
                $doListFile = ".\DoListTask_$($protein)_$number_$i.txt"
                $doListContent = @(
                    "BaseCreate     $i MSAR\$protein`_$number",
                    "BaseAddNewData uncharged\$fileName  Drug",
                    "BaseTraining",
                    "BaseSave",
                    "BaseValidation",
                    "BaseSaveAsSAS",
                    "BaseClose"
                )
                $doListContent | Set-Content -Path $doListFile
                Write-Host "$doListContent"

                # Run the .exe program in a separate, non-blocking process
                try {
                    $process = Start-Process -FilePath $exePath -ArgumentList $doListFile -NoNewWindow -PassThru -RedirectStandardError ".\error.log"
                    $process.WaitForExit()

                    # Append the output to the consolidated output file
                    $output = Get-Content -Path ".\error.log" -ErrorAction SilentlyContinue
                    if ($output) {
                        Add-Content -Path $consolidatedOutputFile -Value "Error running $exePath with $doListFile:`n$output"
                    }
                } catch {
                    Write-Host "Error executing $exePath: $_"
                    Add-Content -Path $consolidatedOutputFile -Value "Error executing $exePath: $_"
                }

                # Optional: Clean up DoListTask file
                Remove-Item -Path $doListFile -Force
            } else {
                Write-Host "Directory $msarDirectory does not exist. Skipping."
            }
        }
    }
}

Write-Host "Task completed. Consolidated output saved to $consolidatedOutputFile."
