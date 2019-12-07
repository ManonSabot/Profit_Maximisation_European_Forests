all: inilialise run_experiments summary_outputs clean

inilialise:
	src/calib_Zroot.sh

run_experiments:
	src/calib_g1.sh
	src/calib_fw.sh

summary_outputs:
	src/postprocessing_scripts/summary_SI.py

clean:
	@echo "Cleaning up all but 10 most recent logs..."
	@ls -tp src/tmp/log.o*| tail -n +11| xargs -r rm --
