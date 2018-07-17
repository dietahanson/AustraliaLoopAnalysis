Running a loop model


There are three main steps:

1. Group species in the network
2. Select which model performs best
3. Run the loop analysis on the best model



*Step 1: Group species in the network

If you have more than 30-40 species in the network, they should be lumped into a smaller number of groups. This can be done by hand (“expertise”/“natural history knowledge”) or using an algorithm such as those in the program n_w. Once you have the species grouped (if needed), format the network as described below.

Output: 

	a) csv table of the grouped network. This should have three columns: 
		1. from: This is the species (or group) that the interaction is coming from. In the case of a trophic interaction, this is the predator.  
		2. to: This is the species (or group) that the interaction is affecting. In the case of a trophic interaction, this is the prey.
		3. type: This is the type of interaction. Current valid types are: predatorprey, habitat, positive, competition.

		Example:

		from,to,type
		bear,rabbit,predatorprey
		rabbit,grass,predatorprey
		rabbit,mouse,competition
		grass,mouse,habitat



*Step 2: Select which model performs best

This step uses the “model_select.R” script to select which model most often gives outcomes of the monitored species (or groups) that match your known outcomes. To run, open the script and change the path for “networktable”, “models”, and “outcomes” to your locations. If you want to change the number of times the models will get compared from the default, change “n.samples”. After running, a plot will get displayed showing the posterior probabilities of model correctness, and a csv table “pp.csv” will get saved with the same probabilities.

Input:

	a) “networktable”: the csv table describing network from step 1

	b) “models”: a csv table describing the different models you want to test. This should have one row for each model and one column for each possible interaction type you want to include. Inclusion of a interaction type in a model is coded with a “1” and exclusion is coded with a “0”. The interaction types in this table must match the interaction types in the network table.

	Example:

	predatorprey, competition, habitat
	1,0,0
	1,0,1
	1,1,1

	c) “outcomes”: a csv table of the known outcomes which will be used to validate the models. No column headers are needed. The first column should contain the species (or groups) for which you have known outcomes, and the second column should contain 1 or -1 to specify the direction of the known outcome. IMPORTANT: the first entry of this table should be the node which is being pressed. So, if the model is being used to predict the effects of removing rabbits from the system, the first line of this file should be “rabbit, -1”. This entry will also be used as one of the validation criteria, that is, if we are testing the effects of removing rabbits from the system, then the model should only be valid if rabbits do indeed decrease.

	Example:

	rabbit, -1
	mouse, -1
	bear, 1

Output:

	a) “pp.csv”: a csv table


*Step 3: Run the loop analysis using the best model

Using the model selected in step 2, run the model using the "loop_model.R" script to run the loop analysis and get predicted outcomes. This takes the same network, outcomes, and press used in step 2 and runs the loop analysis using only the selected model. The inputs are the same "networktable" and "outcomes" as used in step 2. After running, a plot of the predicted outcomes will be displayed, and a csv table "summary" of the same values will get saved.

Input:

	a) “networktable”: the csv table describing network from step 1

	b) “outcomes”: a csv table of the known outcomes used in step 2

Output:

	a) "summary": a csv table showing all the outcomes of the loop analysis--how many times a species/group showed a positive, negative, or neutral trend after the press perturbation.







 
 

