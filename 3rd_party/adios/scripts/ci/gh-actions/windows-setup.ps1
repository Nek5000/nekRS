$ErrorActionPreference = "Stop"

Write-Host "::group::Updating pip"
python3 -m pip install --upgrade pip
Write-Host "::endgroup::"

Write-Host "::group::Installing NumPy"
pip install "numpy>=1.19"
Write-Host "::endgroup::"

Write-Host "::group::Setup CONDA"
$Env:Path += ";$Env:CONDA\condabin"
conda.bat init powershell
conda.bat init bash
Write-Host "::endgroup::"

Write-Host "::group::Installing c-blosc2"
conda.bat install -y conda-forge::c-blosc2
Write-Host "::endgroup::"

if($Env:GH_YML_MATRIX_PARALLEL -eq "ompi")
{
  # This is taken from the MSMPI VCPKG
  $baseurl = "https://download.microsoft.com/download/a/5/2/a5207ca5-1203-491a-8fb8-906fd68ae623"
  $version = "10.1.12498"

  $tempdir    = $Env:RUNNER_TEMP
  $msmpisdk   = Join-Path $tempdir msmpisdk.msi
  $msmpisetup = Join-Path $tempdir msmpisetup.exe

  Write-Host "::group::Downloading Microsoft MPI SDK $version"
  Invoke-WebRequest "$baseurl/msmpisdk.msi" -OutFile $msmpisdk
  Write-Host "::endgroup::"
  Write-Host "::group::Installing Microsoft MPI SDK $version"
  Start-Process msiexec.exe -ArgumentList "/quiet /passive /qn /i $msmpisdk" -Wait
  Write-Host "::endgroup::"

  Write-Host "::group::Downloading Microsoft MPI Runtime $version"

  Invoke-WebRequest "$baseurl/msmpisetup.exe" -OutFile $msmpisetup
  Write-Host "::endgroup::"
  Write-Host "::group::Installing Microsoft MPI Runtime $version"
  Start-Process $msmpisetup -ArgumentList "-unattend" -Wait
  Write-Host "::endgroup::"

  if ($Env:GITHUB_ENV) {
    Write-Host '::group::Adding environment variables to $GITHUB_ENV'
      $envlist = @("MSMPI_BIN", "MSMPI_INC", "MSMPI_LIB32", "MSMPI_LIB64")
      foreach ($name in $envlist) {
        $value = [Environment]::GetEnvironmentVariable($name, "Machine")
          Write-Host "$name=$value"
          Add-Content $Env:GITHUB_ENV "$name=$value"
      }
    Write-Host "::endgroup::"
  }

  if ($Env:GITHUB_PATH) {
    Write-Host '::group::Adding $MSMPI_BIN to $GITHUB_PATH'
      $MSMPI_BIN = [Environment]::GetEnvironmentVariable("MSMPI_BIN", "Machine")
      Add-Content $Env:GITHUB_PATH $MSMPI_BIN
      Write-Host "::endgroup::"
  }

  Write-Host "::group::Installing mpi4py"
  pip install "mpi4py>=1.03"
  Write-Host "::endgroup::"
}

if ($Env:GH_YML_MATRIX_COMPILER -eq "vs2022")
{
  $url = "https://github.com/HDFGroup/hdf5/releases/download/hdf5-1_14_1-2/hdf5-1_14_1-2-win_vs2022.zip"

  $tempdir    = $Env:RUNNER_TEMP
  $hdf5Zip    = Join-Path $tempdir hdf5-1_14_1-2-win_vs2022.zip
  $installDir  = "C:\hdf5"
  $extractDir = Join-Path $tempdir hdf5
  $innerZip   = Join-Path $extractDir HDF5-1.14.2.1-win64.zip

  Write-Host "::group::Downloading windows HDF5 release $url"
  Invoke-WebRequest $url -OutFile $hdf5Zip
  Write-Host "::endgroup::"

  Write-Host "::group::Extracting outer windows HDF5 release $hdf5Zip"
  Expand-Archive $hdf5Zip -DestinationPath $tempDir -Force
  Write-Host "::endgroup::"

  Write-Host "::group::Extracting inner windows HDF5 release $innerZip"
  Expand-Archive $innerZip -DestinationPath $installDir -Force
  Write-Host "::endgroup::"

  # Now the HDF5_DIR can be set to: C:\hdf5\HDF5-1.14.2.1-win64\cmake
}
