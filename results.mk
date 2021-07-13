all: run_experiments plot_results

run_experiments:
	src/predict_kmax.sh -A
	src/var_solver_accuracy.sh

plot_results:
	src/Sabot_et_al_2019_all_plots.sh

clean:
	@echo "Cleaning up all but 10 most recent logs..."
	@ls -tp src/tmp/log.o*| tail -n +11| xargs -r rm --
