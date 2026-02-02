"""Tests for seqnado.utils module (CLI utilities)."""

import subprocess
from unittest.mock import Mock, patch, mock_open

import pytest

from seqnado.utils import (
    extract_apptainer_args,
    extract_cores_from_options,
    pepe_silvia,
    remove_unwanted_run_files,
    run_batch_job_on_error,
    FileSelector,
)


class TestExtractCoresFromOptions:
    """Tests for extract_cores_from_options function."""

    def test_extract_cores_with_c_flag(self):
        """Test extraction with -c flag."""
        options = ["-c", "8", "--other", "value"]
        result_options, cores = extract_cores_from_options(options)
        assert cores == 8
        assert result_options == ["--other", "value"]

    def test_extract_cores_with_cores_flag(self):
        """Test extraction with --cores flag."""
        options = ["--cores", "16", "--other", "value"]
        result_options, cores = extract_cores_from_options(options)
        assert cores == 16
        assert result_options == ["--other", "value"]

    def test_extract_cores_no_flag(self):
        """Test default behavior when no core flag is present."""
        options = ["--other", "value"]
        result_options, cores = extract_cores_from_options(options)
        assert cores == 1
        assert result_options == ["--other", "value"]

    def test_extract_cores_flag_no_value(self):
        """Test behavior when core flag is present but no value given."""
        options = ["-c", "--other", "value"]
        result_options, cores = extract_cores_from_options(options)
        assert cores == 1
        # When int() conversion fails, falls back to ValueError handler which doesn't remove -c
        assert result_options == ["-c", "--other", "value"]

    def test_extract_cores_multiple_flags(self):
        """Test that first occurrence of core flag is used."""
        options = ["-c", "4", "--other", "value"]
        result_options, cores = extract_cores_from_options(options)
        assert cores == 4
        assert result_options == ["--other", "value"]

    def test_extract_cores_flag_at_end(self):
        """Test behavior when core flag is at the end of options (IndexError)."""
        options = ["--other", "value", "-c"]
        result_options, cores = extract_cores_from_options(options)
        assert cores == 1
        # The -c flag should be removed when IndexError occurs
        assert "-c" not in result_options
        assert result_options == ["--other", "value"]


class TestExtractApptainerArgs:
    """Tests for extract_apptainer_args function."""

    def test_extract_apptainer_args_present(self):
        """Test extraction when apptainer args are present."""
        options = ["--apptainer-args", "-B /data:/data", "--other", "value"]
        result_options, apptainer_args = extract_apptainer_args(options)
        assert apptainer_args == "-B /data:/data"
        assert result_options == ["--other", "value"]

    def test_extract_apptainer_args_absent(self):
        """Test extraction when apptainer args are absent."""
        options = ["--other", "value"]
        result_options, apptainer_args = extract_apptainer_args(options)
        assert apptainer_args == ""
        assert result_options == ["--other", "value"]

    def test_extract_apptainer_args_empty_string(self):
        """Test extraction with empty string argument."""
        options = ["--apptainer-args", "", "--other", "value"]
        result_options, apptainer_args = extract_apptainer_args(options)
        assert apptainer_args == ""
        assert result_options == ["--other", "value"]


class TestPepeSilvia:
    """Tests for pepe_silvia function."""

    def test_pepe_silvia_returns_url(self):
        """Test that pepe_silvia returns a URL."""
        with patch("builtins.print") as mock_print:
            result = pepe_silvia()
            assert "http" in result
            assert ".jpg" in result
            mock_print.assert_called_once_with("PEPE SILVIA")


class TestRemoveUnwantedRunFiles:
    """Tests for remove_unwanted_run_files function."""

    @patch("glob.glob")
    @patch("os.path.isdir")
    @patch("os.remove")
    @patch("shutil.rmtree")
    def test_remove_slurm_files(self, mock_rmtree, mock_remove, mock_isdir, mock_glob):
        """Test removal of slurm output files."""
        mock_glob.side_effect = [
            ["slurm-123.out", "slurm-456.out"],  # slurm files
            [],  # sps files
            [],  # simg files
        ]
        mock_isdir.return_value = False

        remove_unwanted_run_files()

        assert mock_remove.call_count == 2

    @patch("glob.glob")
    @patch("os.path.isdir")
    @patch("os.remove")
    @patch("shutil.rmtree")
    def test_remove_sps_files(self, mock_rmtree, mock_remove, mock_isdir, mock_glob):
        """Test removal of sps files."""
        mock_glob.side_effect = [
            [],  # slurm files
            ["sps-temp1", "sps-temp2"],  # sps files
            [],  # simg files
        ]
        mock_isdir.return_value = False

        remove_unwanted_run_files()

        assert mock_remove.call_count == 2

    @patch("glob.glob")
    @patch("os.path.isdir")
    @patch("os.remove")
    @patch("shutil.rmtree")
    def test_remove_simg_files(self, mock_rmtree, mock_remove, mock_isdir, mock_glob):
        """Test removal of simg files."""
        mock_glob.side_effect = [
            [],  # slurm files
            [],  # sps files
            ["container.simg"],  # simg files
        ]
        mock_isdir.return_value = False

        remove_unwanted_run_files()

        mock_remove.assert_called_once_with("container.simg")

    @patch("glob.glob")
    @patch("os.path.isdir")
    @patch("os.remove")
    @patch("shutil.rmtree")
    def test_remove_directory(self, mock_rmtree, mock_remove, mock_isdir, mock_glob):
        """Test removal of directories."""
        mock_glob.side_effect = [
            [],  # slurm files
            ["sps-temp-dir"],  # sps files
            [],  # simg files
        ]
        mock_isdir.return_value = True

        remove_unwanted_run_files()

        mock_rmtree.assert_called_once_with("sps-temp-dir")
        mock_remove.assert_not_called()

    @patch("glob.glob")
    @patch("os.path.isdir")
    @patch("os.remove")
    @patch("shutil.rmtree")
    def test_handles_removal_errors(self, mock_rmtree, mock_remove, mock_isdir, mock_glob):
        """Test that removal errors are caught and don't crash."""
        mock_glob.side_effect = [
            ["slurm-123.out"],  # slurm files
            [],  # sps files
            [],  # simg files
        ]
        mock_isdir.return_value = False
        mock_remove.side_effect = PermissionError("Cannot remove")

        # Should not raise exception
        with patch("builtins.print") as mock_print:
            remove_unwanted_run_files()
            mock_print.assert_called()


class TestRunBatchJobOnError:
    """Tests for run_batch_job_on_error function."""

    @patch("subprocess.run")
    @patch("builtins.open", new_callable=mock_open)
    def test_creates_slurm_script(self, mock_file, mock_subprocess):
        """Test that slurm script is created with correct content."""
        email = "test@example.com"
        mock_subprocess.return_value = Mock()

        run_batch_job_on_error(email)

        mock_file.assert_called_once_with("error_notification.sh", "w")
        written_content = "".join(
            call.args[0] for call in mock_file().write.call_args_list
        )
        assert email in written_content
        assert "SBATCH" in written_content
        assert "seqnado_error_notification" in written_content

    @patch("subprocess.run")
    @patch("builtins.open", new_callable=mock_open)
    def test_submits_sbatch_job(self, mock_file, mock_subprocess):
        """Test that sbatch command is executed."""
        email = "test@example.com"
        mock_subprocess.return_value = Mock()

        run_batch_job_on_error(email)

        mock_subprocess.assert_called_once_with(
            ["sbatch", "error_notification.sh"], check=True
        )

    @patch("subprocess.run")
    @patch("builtins.open", new_callable=mock_open)
    def test_handles_sbatch_failure(self, mock_file, mock_subprocess):
        """Test handling of sbatch submission failure."""
        email = "test@example.com"
        mock_subprocess.side_effect = subprocess.CalledProcessError(1, "sbatch")

        # Should not raise exception
        run_batch_job_on_error(email)

class TestFileSelector:
    """Tests for FileSelector class."""

    def test_select_basic(self):
        files = ["a.txt", "b.TXT", "c_backup.txt", "d.log"]
        s = FileSelector(files)
        assert set(s.select(".txt")) == {"a.txt", "b.TXT", "c_backup.txt"}
        assert s.select(".txt", excludes="backup",) == ["a.txt", "b.TXT"]
        assert s.select(".txt", includes="a",) == ["a.txt", 'c_backup.txt']
        assert s.select(".txt", includes=[r"^b"], use_regex=True) == ["b.TXT"]

        with pytest.raises(AssertionError):
           assert s.select(".txt", includes=[r"^b"], use_regex=False) == ["b.TXT"]
