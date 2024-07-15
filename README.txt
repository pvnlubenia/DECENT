========================================================

  DECENT: DEcompositions of Chemical rEaction NeTworks

========================================================

MATLAB was used to develop the functions used here.

This package contains four functions:

   1. FIID_check
   2. FID_check
   3. FIID
   4. FID

All of these deal with finest decompositions and the network numbers of the resulting subnetworks. Use the functions with "_check" to instantly generate values for only three network numbers (see below). Use the "pure" functions to generate complete tables of network numbers.

The functions FIID_check and FID_check return the nontrivial finest incidence independent decomposition and nontrivial finest independent decomposition, respectively, of a chemical reaction network (CRN), if it exists. If no such decomposition exists, a message appears saying so. They also show a table summarizing the deficiency, reactant deficiency, and reactant rank values of the entire network and of the subnetworks. The sum of the subnetworks' network numbers, which may provide additional insights, also appear at the end of the table.

Similar outputs appear for the functions FIID and FID, but a complete table summarizing all the network numbers appear.
For FIID_check and FIID, the output variables 'model', 'I_a', 'G', 'P', and 'results' allow the user to view the following, respectively:

   - Complete network with all the species listed in the 'species' field of the structure 'model'
   - Incidence matrix of the network
   - Undirected graph of I_a
   - Partitions representing the decomposition of the reactions
   - Complete table of network numbers

For FID_check and FID, 'I_a' is replaced by 'R', the matrix of reaction vectors of the network (consequently, 'G' becomes the undirected graph of R).

The algorithm for the decomposition in FID_check and FID are based on [3]. A modification of the algorithm to apply to FIID_check and FID were based on [2]. Concepts in the network numbers are based on [1] while their codes are largely derived from [7].



=================================
How to fill out 'model' structure
=================================

'model' is the input for the functions in this package. It is a structure, representing the CRN, with the following fields (the kinetics of the network is not needed):

   - id: name of the model
   - species: a list of all species in the network; this is left blank since incorporated into the function is a step which compiles all species used in the model
   - reaction: a list of all reactions in the network, each with the following subfields:
        - id: a string representing the reaction
        - reactant: has the following further subfields:
             - species: a list of strings representing the species in the reactant complex
             - stoichiometry: a list of numbers representing the stoichiometric coefficient of each species in the reactant complex (listed in the same order of the species)
        - product: has the following further subfields:
             - species: a list of strings representing the species in the product complex
             - stoichiometry: a list of numbers representing the stoichiometric coefficient of each species in the product complex (listed in the same order of the species)
        - reversible: has the value true or false indicating if the reaction is reversible or not, respectively

To fill out the 'model' structure, write a string for 'model.id': this is just to put a name to the network. To add the reactions to the network, use the function addReaction where the output is 'model'. addReaction is developed to make the input of reactions of the CRN easier than the input in [7]:

   addReaction
      - OUTPUT: Returns a structure called 'model' with added field 'reaction' with subfields 'id', 'reactant', 'product', 'reversible', and 'kinetic'. The output variable 'model' allows the user to view the network with the added reaction.
      - INPUTS
           - model: a structure, representing the CRN
           - id: visual representation of the reaction, e.g., reactant -> product (string)
           - reactant_species: species of the reactant complex (cell)
           - reactant_stoichiometry: stoichiometry of the species of the reactant complex (cell)
           - reactant_kinetic: kinetic orders of the species of the reactant complex (array)
           - product_species: species of the product complex (cell)
           - product_stoichiometry: stoichiometry of the species of the product complex (cell)
           - product_kinetic: "kinetic orders" of the species of the product complex, if the reaction is reversible (array); if the reaction in NOT reversible, leave blank
           - reversible: logical; whether the reaction is reversible or not (true or false)
      * Make sure the function addReaction is in the same folder/path being used as the current working directory.



========
Examples
========

8 examples are included in this folder:

   - Example 1: Baccam influenza virus model [3]

   - Example 2: Handel influenza virus model [3]

   - Example 3: Generalized mass action model of purine metabolism in man [3]

   - Example 4: Insulin metabolic signaling network [5]

   - Example 5: Insulin resistance model [6]

   - Example 6: Schmitz Wnt signaling network [4]

   - Example 7: Feinberg augmented Lee Wnt signaling network [4]

   - Example 8: MacLean Wnt signaling network [4]



===================
Contact Information
===================

For questions, comments, and suggestions, feel free to contact me at pvnlubenia@yahoo.co.uk.


- Patrick Lubenia (15 July 2024)



==========
References
==========

   [1] Arceo C, Jose E, Lao A, Mendoza E (2017) Reactant subspaces and kinetics of chemical reaction networks. J Math Chem 56(5):395-422. https://doi.org/10.1007/s10910-017-0809-x

   [2] Hernandez B, Amistas D, De la Cruz R, Fontanil L, de los Reyes V A, Mendoza E (2022) Independent, incidence independent and weakly reversible decompositions of chemical reaction networks. MATCH Commun Math Comput Chem, 87(2):367--396. https://doi.org/10.46793/match.87-2.367H

   [3] Hernandez B, De la Cruz R (2021) Independent decompositions of chemical reaction networks. Bull Math Biol 83(76):1Ð23. https://doi.org/10.1007/s11538-021-00906-3

   [4] Hernandez B, Lubenia P, Mendoza E (2024) Embedding-based comparison of reaction networks of Wnt signaling. MATCH Commun Math Comput Chem, 93:223-245. https://doi.org/10.46793/match.93-1.223H

   [5] Lubenia P, Mendoza E, Lao A (2022) Reaction network analysis of metabolic insulin signaling. Bull Math Biol, 84(129):1Ð12. https://doi.org/10.1007/s11538-022-01087-3

   [6] Lubenia P, Mendoza E, Lao A (2024) Comparative analysis of kinetic realizations of insulin signaling. J Theor Biol, 577:1Ð12. https://doi.org/10.1016/j.jtbi.2023.111672

   [7] Soranzo N, Altafini C (2009) ERNEST: a toolbox for chemical reaction network theory. Bioinform 25(21):2853-2854. https://doi.org/10.1093/bioinformatics/btp513