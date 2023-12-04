$processorCount = (Get-WmiObject -Class Win32_Processor | Measure-Object -Property NumberOfLogicalProcessors -Sum).Sum

# Define the script block that will be executed in each job
$MINI=1
$MAXI=3
$MINJ=73
$MAXJ=400

$scriptBlock = {
    param($i, $j)
    Set-Location "E:\GitHub\GGM\Gof"
$logFile = "E:\Programming\LocalExperiments\GGM\GoF\MyLogs\log_i_${i}_j_${j}.txt" 
    # Run the R script with arguments
    Rscript "Experiment.R" $i $j > $logFile
}

$jobs = @()

# Loop through each combination of i and j
for ($i = $MINI; $i -le $MAXI; $i++) {
    for ($j = $MINJ; $j -le $MAXJ; $j++) {
	
	while (($jobs | Where-Object { $_.State -eq 'Running' } | Measure-Object).Count -ge $processorCount) {
	        Start-Sleep -Seconds 1
	    }
	
        # Start a job with the current values of i and j
         $jobs += Start-Job -ScriptBlock $scriptBlock -ArgumentList $i, $j
    }
}

# Wait for all jobs to complete
$jobs | Wait-Job

# Retrieve and display the results of all jobs
$jobs | Receive-Job

# Remove all completed jobs
$jobs | Remove-Job
