#! /usr/bin/env python

import os
import matplotlib as mpl
if os.environ.get('DISPLAY', '') == '':
    mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pyemma
from pyemma.util.contexts import settings
import mdtraj as md
from argparse import ArgumentParser
import mdtraj


# Main - used to load data and execute functions
def main():
    # Get commandline parser
    args = cmdlineparse()

    # Parse residue pairs for sidechain distance calculation
    residue_pairs_sc = []
    if args.residue_pair_sc_files:
        for file in args.residue_pair_sc_files:
            with open(file, 'r') as infile:
                file_content = infile.readlines()
                for line in file_content:
                    residue_pairs_sc.append([int(x) for x in line.strip().split(',')])
    # Parse residue pairs
    residue_pairs_ch = []
    if args.residue_pair_closest_files:
        for file in args.residue_pair_closest_files:
            with open(file, 'r') as infile:
                file_content = infile.readlines()
                for line in file_content:
                    residue_pairs_ch.append([int(x) for x in line.strip().split(',')])

    # Parse atom pairs
    atom_pairs = []
    if args.atom_pairs_files:
        for file in args.atom_pairs_files:
            with open(file, 'r') as infile:
                file_content = infile.readlines()
                for line in file_content:
                    atom_pairs.append([int(x) for x in line.strip().split(',')])

    # Parse atom iterables A
    atom_iterables_a = []
    if args.atom_iterables_a_files:
        for file in args.atom_iterables_a_files:
            with open(file, 'r') as infile:
                file_content = infile.readlines()
                for line in file_content:
                    atom_iterables_a.append([int(x) for x in line.strip().split(',')[:-1]])

    # Parse atom iterables B
    atom_iterables_b = []
    if args.atom_iterables_b_files:
        for file in args.atom_iterables_b_files:
            with open(file, 'r') as infile:
                file_content = infile.readlines()
                for line in file_content:
                    print(line)
                    atom_iterables_b.append([int(x) for x in line.strip().split(',')[:-1]])

    # Parse atom contacts
    atom_contacts = []
    if args.atom_contacts_files:
        for file in args.atom_contacts_files:
            with open(file, 'r') as infile:
                file_content = infile.readlines()
                for line in file_content:
                    atom_contacts.append([int(x) for x in line.strip().split(',')])

    # Parse dihedral atoms
    dihedrals = []
    if args.dihedral_atoms_files:
        for file in args.dihedral_atoms_files:
            with open(file, 'r') as infile:
                file_content = infile.readlines()
                for line in file_content:
                    dihedrals.append([int(x) for x in line.strip().split(',')])

    # Parse atom groups
    atom_groups = []
    if args.atom_groups_files:
        for file in args.atom_groups_files:
            with open(file, 'r') as infile:
                file_content = infile.readlines()
                for line in file_content:
                    atom_groups.append([int(x) for x in line.replace("[", "").replace("]", "").strip().split(", ")])

    # Parse atom groups
    group_pairs = []
    if args.group_pairs_file:
        with open(args.group_pairs_file, 'r') as infile:
            file_content = infile.readlines()
            for line in file_content:
                group_pairs.append([int(x) for x in line.strip().split(',')])

    # Run Featurize
    features = Featurize(args, residue_pairs_sc, residue_pairs_ch, atom_pairs, atom_iterables_a, atom_iterables_b,
                         atom_contacts, dihedrals, atom_groups, group_pairs)

    # Store features to reader
    if args.load_tica_for_projection or args.calc_tica_for_projection:
        reader = pyemma.coordinates.load(args.trajfiles, features=features)  # transform only works on loaded data
    else:
        reader = pyemma.coordinates.source(args.trajfiles, features=features)

    # Load a TICA trajectory to project new input data onto
    if args.load_tica_for_projection:
        tica_ref = pyemma.load(str(args.load_tica_for_projection[0]), model_name=str(args.checkpoint_load[1]) + "_tica")

    # Compute a reference TICA trajectory to project new data onto
    elif args.calc_tica_for_projection:
        ref_features = FeaturizeRef(args)
        ref_reader = pyemma.coordinates.load(args.calc_tica_for_projection[1],
                                           features=ref_features)  # transform only works on loaded data
        # ref_reader = pyemma.coordinates.source(args.calc_tica_for_projection[1], features=ref_features)
        tica_ref = pyemma.coordinates.tica(ref_reader, lag=int(args.lag[0]), var_cutoff=float(args.tica),
                                           stride=int(args.dim_red_stride))

    # Load pre-saved features
    if args.checkpoint_load:
        feat_traj = reader.get_output(stride=int(args.stride))
        if args.tica:
            tica = pyemma.load(str(args.checkpoint_load[0]), model_name=str(args.checkpoint_load[1]) + "_tica")
            if args.kmeans:
                cluster = pyemma.coordinates.cluster_kmeans(tica, k=int(args.k_clusters), max_iter=int(args.k_iter),
                                                            stride=int(args.dim_red_stride))
            elif args.regspace:
                cluster = pyemma.coordinates.cluster_regspace(tica, dmin=float(args.regspace_dmin),
                                                              stride=int(args.dim_red_stride))
        if args.pca:
            pca = pyemma.load(str(args.checkpoint_load[0]), model_name=str(args.checkpoint_load[1]) + "_pca")
            if args.kmeans:
                cluster = pyemma.coordinates.cluster_kmeans(pca, k=int(args.k_clusters), max_iter=int(args.k_iter),
                                                            stride=int(args.dim_red_stride))
            elif args.regspace:
                cluster = pyemma.coordinates.cluster_regspace(pca, dmin=float(args.regspace_dmin),
                                                              stride=int(args.dim_red_stride))
        msm = pyemma.load(str(args.checkpoint_load[0]), model_name=str(args.checkpoint_load[1]) + "_msm")
    #         cluster = pyemma.load(str(args.checkpoint_load[0]), model_name=str(args.checkpoint_load[1])+"_cluster")

    elif args.load_tica_for_projection or args.calc_tica_for_projection:
        tica = tica_ref.transform(reader)
        feat_traj = tica
        if args.kmeans:
            cluster = pyemma.coordinates.cluster_kmeans(tica, k=int(args.k_clusters), max_iter=int(args.k_iter),
                                                        stride=int(args.dim_red_stride))
        elif args.regspace:
            cluster = pyemma.coordinates.cluster_regspace(tica, dmin=float(args.regspace_dmin),
                                                          stride=int(args.dim_red_stride))

        # Build a Markov model
        if args.kmeans or args.regspace:
            if args.hidden:
                nstates = int(args.n_metastable_states)
                msm = pyemma.msm.estimate_hidden_markov_model(cluster.dtrajs, nstates, lag=int(args.lag[0]),
                                                              dt_traj=str(args.dt))
            else:
                msm = pyemma.msm.estimate_markov_model(cluster.dtrajs, lag=int(args.lag[0]), dt_traj=str(args.dt))
        else:
            if args.hidden:
                nstates = int(args.n_metastable_states)
                msm = pyemma.msm.estimate_hidden_markov_model(feat_traj, nstates, lag=int(args.lag[0]),
                                                              dt_traj=str(args.dt))
            else:
                msm = pyemma.msm.estimate_markov_model(feat_traj, lag=int(args.lag[0]), dt_traj=str(args.dt))

    # Perform new analysis
    else:
        # Generate the feature trajectory
        if args.tica:
            print("Performing new TICA analysis with a variance cutoff of " + str(args.tica))
            tica = pyemma.coordinates.tica(reader, lag=int(args.lag[0]), var_cutoff=float(args.tica),
                                           stride=int(args.dim_red_stride))
            feat_traj = tica.get_output(stride=int(args.stride))
            if args.checkpoint_save:
                tica.save(str(args.checkpoint_save[0]), model_name=str(args.checkpoint_save[1]) + "_tica",
                          overwrite=True)
            if args.kmeans:
                cluster = pyemma.coordinates.cluster_kmeans(tica, k=int(args.k_clusters), max_iter=int(args.k_iter),
                                                            stride=int(args.dim_red_stride))
            elif args.regspace:
                cluster = pyemma.coordinates.cluster_regspace(tica, dmin=float(args.regspace_dmin),
                                                              stride=int(args.dim_red_stride))
        elif args.pca:
            pca = pyemma.coordinates.pca(reader, var_cutoff=float(args.pca), stride=int(args.dim_red_stride))
            feat_traj = pca.get_output(stride=int(args.stride))
            if args.checkpoint_save:
                pca.save(str(args.checkpoint_save[0]), model_name=str(args.checkpoint_save[1]) + "_pca", overwrite=True)
            if args.kmeans:
                cluster = pyemma.coordinates.cluster_kmeans(pca, k=int(args.k_clusters), max_iter=int(args.k_iter),
                                                            stride=int(args.dim_red_stride))
            elif args.regspace:
                cluster = pyemma.coordinates.cluster_regspace(pca, dmin=float(args.regspace_dmin),
                                                              stride=int(args.dim_red_stride))
        elif args.kmeans or args.regspace:
            if args.kmeans:
                cluster = pyemma.coordinates.cluster_kmeans(reader.get_output(stride=int(args.dim_red_stride)),
                                                            k=int(args.k_clusters), max_iter=int(args.k_iter))
            elif args.regspace:
                cluster = pyemma.coordinates.cluster_regspace(reader.get_output(stride=int(args.dim_red_stride)),
                                                              dmin=float(args.regspace_dmin))
            feat_traj = reader.get_output(stride=int(args.stride))
        else:
            feat_traj = reader.get_output(stride=int(args.stride))

        # Save clusters if made
        #         if (args.kmeans or args.regspace) and args.checkpoint_save:
        #             cluster.save(str(args.checkpoint_save[0]), model_name=str(args.checkpoint_save[1])+"_cluster", overwrite=True, save_streaming_chain=False) # If I turn on, saving fails, but if I turn off, Plot fails
        #
        # Build a Markov model
        if args.kmeans or args.regspace:
            if args.hidden:
                nstates = int(args.n_metastable_states)
                msm = pyemma.msm.estimate_hidden_markov_model(cluster.dtrajs, nstates, lag=int(args.lag[0]),
                                                              dt_traj=str(args.dt))
            else:
                msm = pyemma.msm.estimate_markov_model(cluster.dtrajs, lag=int(args.lag[0]), dt_traj=str(args.dt))
        else:
            if args.hidden:
                nstates = int(args.n_metastable_states)
                msm = pyemma.msm.estimate_hidden_markov_model(feat_traj, nstates, lag=int(args.lag[0]),
                                                              dt_traj=str(args.dt))
            else:
                msm = pyemma.msm.estimate_markov_model(feat_traj, lag=int(args.lag[0]), dt_traj=str(args.dt))

        # Check reversible connectivity
        if not args.hidden:
            print('fraction of states used = {:f}'.format(msm.active_state_fraction))
            print('fraction of counts used = {:f}'.format(msm.active_count_fraction))

        # Perform CK test with a given lag time
        if args.ck_test:
            pyemma.plots.plot_cktest(msm.cktest(int(args.ck_test)), units=str(args.dt[1]));
            plt.savefig(str(args.output_prefix) + ".CK_Test.png", dpi=300)
            plt.clf()

        # Save our model and everything that went into it
        if args.checkpoint_save:
            msm.save(str(args.checkpoint_save[0]), model_name=str(args.checkpoint_save[1]) + "_msm", overwrite=True)

    # Compute implied time scales
    if len(args.lag) > 1 and not (args.load_tica_for_projection or args.calc_tica_for_projection):
        if args.kmeans or args.regspace:
            if args.kmeans:
                cluster = pyemma.coordinates.cluster_kmeans(reader.get_output(stride=int(args.stride)),
                                                            k=int(args.k_clusters), max_iter=int(args.k_iter),
                                                            stride=int(args.dim_red_stride))
            elif args.regspace:
                cluster = pyemma.coordinates.cluster_regspace(reader.get_output(stride=int(args.stride)),
                                                              dmin=float(args.regspace_dmin),
                                                              stride=int(args.dim_red_stride))
            its = pyemma.msm.its(cluster.dtrajs, lags=[int(x) for x in args.lag], nits=int(args.nits), errors='bayes')
        else:
            its = pyemma.msm.its(reader.get_output(stride=int(args.stride)), lags=[int(x) for x in args.lag],
                                 nits=int(args.nits), errors='bayes')
        pyemma.plots.plot_implied_timescales(its, ylog=True);
        plt.savefig(str(args.output_prefix) + ".ITS.png", dpi=300)
        plt.clf()
        # don't do additional calculations until ITS is viewed
        return 0
    else:
        its = []

    # Run Plot
    print("Sample microstates from each macrostate!")
    SampleStates(msm, args)
    print("Plot free energy landscapes!")
    Plot(feat_traj, features, cluster, its, msm, args)


# Generate feature trajectories
def Featurize(ARGS, RES_PAIRS_SC, RES_PAIRS_CH, ATOM_PAIRS, ATOM_ITERABLES_A, ATOM_ITERABLES_B, ATOM_CONTACTS,
              DIHEDRALS, ATOM_GROUPS, GROUP_PAIRS):
    # Load topology file and generate feature object
    features = pyemma.coordinates.featurizer(ARGS.topfile)

    # Add residue pair sidechain features
    for pair in range(0, len(RES_PAIRS_SC)):
        features.add_residue_mindist([RES_PAIRS_SC[pair]], scheme='sidechain-heavy', threshold=None,
                                     periodic=ARGS.periodic)

    # Add residue pair closest heavy atom features
    for pair in range(0, len(RES_PAIRS_CH)):
        features.add_residue_mindist([RES_PAIRS_CH[pair]], scheme='closest-heavy', threshold=None,
                                     periodic=ARGS.periodic)

    # Add dihedral features
    for dihedrals in range(0, len(DIHEDRALS)):
        features.add_dihedrals([DIHEDRALS[dihedrals]], deg=True, cossin=False, periodic=ARGS.periodic)

    # Add atom pair features
    for pair in range(0, len(ATOM_PAIRS)):
        features.add_distances([ATOM_PAIRS[pair]], periodic=ARGS.periodic)

    # Add atom iterables distance features
    if ARGS.atom_iterables_a_files and ARGS.atom_iterables_b_files:
        features.add_distances(ATOM_ITERABLES_A[0], indices2=ATOM_ITERABLES_B[0], periodic=ARGS.periodic)
    elif ARGS.atom_iterables_a_files:
        features.add_distances(ATOM_ITERABLES_A[0], periodic=ARGS.periodic)

    # Add atom contacts
    for pair in range(0, len(ATOM_CONTACTS)):
        features.add_contacts([ATOM_CONTACTS[pair]], threshold=0.35, periodic=ARGS.periodic)

    # Add geometric center atom group distances
    if ARGS.group_pairs_file:
        GROUP_PAIRS = np.array(GROUP_PAIRS)  # must be an nx2 ndarray
        features.add_group_mindist(ATOM_GROUPS, group_pairs=GROUP_PAIRS, threshold=None, periodic=ARGS.periodic)

    # Add dihedral features
    if ARGS.backbone_torsions:
        features.add_backbone_torsions(selstr=str(ARGS.backbone_torsions), deg=True, cossin=True,
                                       periodic=ARGS.periodic)

    # Add sidechain features
    if ARGS.sidechain_torsions:
        features.add_sidechain_torsions(selstr=str(ARGS.backbone_torsions), deg=True, cossin=True, which='all',
                                        periodic=ARGS.periodic)

    # Print shit and return
    print("Input features: \n")
    print(features.describe())
    print("End input features \n")
    return features


# Generate feature trajectories
def FeaturizeRef(ARGS):
    # Load topology file and generate feature object
    RES_PAIRS_SC, RES_PAIRS_CH, ATOM_PAIRS, DIHEDRALS, ATOM_GROUPS, GROUP_PAIRS = LoadRefData(ARGS)
    features = pyemma.coordinates.featurizer(ARGS.calc_tica_for_projection[0])  # ref topfile

    # Add residue pair sidechain features
    for pair in range(0, len(RES_PAIRS_SC)):
        features.add_residue_mindist([RES_PAIRS_SC[pair]], scheme='sidechain-heavy', threshold=None,
                                     periodic=ARGS.periodic)

    # Add residue pair closest heavy atom features
    for pair in range(0, len(RES_PAIRS_CH)):
        features.add_residue_mindist([RES_PAIRS_CH[pair]], scheme='closest-heavy', threshold=None,
                                     periodic=ARGS.periodic)

    # Add dihedral features
    for dihedrals in range(0, len(DIHEDRALS)):
        features.add_dihedrals([DIHEDRALS[dihedrals]], deg=True, cossin=False, periodic=ARGS.periodic)

    # Add atom pair features
    for pair in range(0, len(ATOM_PAIRS)):
        features.add_distances([ATOM_PAIRS[pair]], periodic=ARGS.periodic)

    # Add geometric center atom group distances
    if ARGS.group_pairs_file:
        GROUP_PAIRS = np.array(GROUP_PAIRS)  # must be an nx2 ndarray
        features.add_group_mindist(ATOM_GROUPS, group_pairs=GROUP_PAIRS, threshold=None, periodic=ARGS.periodic)

    # Add dihedral features
    if ARGS.backbone_torsions:
        features.add_backbone_torsions(deg=True, cossin=True, periodic=ARGS.periodic)

    # Add sidechain features
    if ARGS.sidechain_torsions:
        features.add_sidechain_torsions(deg=True, cossin=True, which='all', periodic=ARGS.periodic)

    # Print shit and return
    print("Reference features: \n")
    print(features.describe())
    print("End reference features \n")
    return features


def LoadRefData(args):
    # Parse residue pairs for sidechain distance calculation
    ref_residue_pairs_sc = []
    if args.ref_residue_pair_sc_files:
        for file in args.ref_residue_pair_sc_files:
            with open(file, 'r') as infile:
                file_content = infile.readlines()
                for line in file_content:
                    ref_residue_pairs_sc.append([int(x) for x in line.strip().split(',')])
    # Parse residue pairs
    ref_residue_pairs_ch = []
    if args.ref_residue_pair_closest_files:
        for file in args.ref_residue_pair_closest_files:
            with open(file, 'r') as infile:
                file_content = infile.readlines()
                for line in file_content:
                    ref_residue_pairs_ch.append([int(x) for x in line.strip().split(',')])

    # Parse dihedral atoms
    ref_dihedrals = []
    if args.ref_dihedral_atoms_files:
        for file in args.ref_dihedral_atoms_files:
            with open(file, 'r') as infile:
                file_content = infile.readlines()
                for line in file_content:
                    ref_dihedrals.append([int(x) for x in line.strip().split(',')])

    # Parse atom pairs
    ref_atom_pairs = []
    if args.ref_atom_pairs_files:
        for file in args.ref_atom_pairs_files:
            with open(file, 'r') as infile:
                file_content = infile.readlines()
                for line in file_content:
                    ref_atom_pairs.append([int(x) for x in line.strip().split(',')])

    # Parse atom groups
    ref_atom_groups = []
    if args.ref_atom_groups_files:
        for file in args.ref_atom_groups_files:
            with open(file, 'r') as infile:
                file_content = infile.readlines()
                for line in file_content:
                    ref_atom_groups.append([int(x) for x in line.replace("[", "").replace("]", "").strip().split(", ")])

    # Parse atom groups
    ref_group_pairs = []
    if args.ref_group_pairs_file:
        with open(args.ref_group_pairs_file, 'r') as infile:
            file_content = infile.readlines()
            for line in file_content:
                ref_group_pairs.append([int(x) for x in line.strip().split(',')])

    # return ref data
    return ref_residue_pairs_sc, ref_residue_pairs_ch, ref_atom_pairs, ref_dihedrals, ref_atom_groups, ref_group_pairs


def SampleStates(MSM, ARGS):
    "Sample from metastable states"
    nstates = int(ARGS.n_metastable_states)
    if not ARGS.hidden:
        MSM.pcca(nstates)

    # Print stationary distributions
    if ARGS.hidden:
        print(MSM.pi)
        ofile = open(str(ARGS.output_prefix) + ".stationary_distributions.dat", "w+")
        for i, s in enumerate(MSM.metastable_sets):
            print(MSM.pi[i])
            ofile.write('π_{} = {:f}'.format(i, MSM.pi[i].sum()))
            ofile.write('\n')
        ofile.close()
    else:
        ofile = open(str(ARGS.output_prefix) + ".stationary_distributions.dat", "w+")
        for i, s in enumerate(MSM.metastable_sets):
            ofile.write('π_{} = {:f}'.format(i, MSM.pi[s].sum()))
            ofile.write('\n')
        ofile.close()

    cmap = mpl.cm.get_cmap('viridis', nstates)
    index = int(0)
    if ARGS.hidden:
        for idist in MSM.sample_by_observation_probabilities(int(ARGS.n_samples_each_state)):
            my_samples = pyemma.coordinates.save_traj(ARGS.trajfiles, idist,
                                                      str(ARGS.output_prefix) + "." + str(index) + ".metastable.nc",
                                                      top=str(ARGS.topfile), image_molecules=bool(ARGS.periodic))
            index = index + 1
        if ARGS.ngl:
            VisualizeMetastable(my_samples, cmap, selection='backbone')
    else:
        for idist in MSM.sample_by_distributions(MSM.metastable_distributions, int(ARGS.n_samples_each_state)):
            my_samples = pyemma.coordinates.save_traj(ARGS.trajfiles, idist,
                                                      str(ARGS.output_prefix) + "." + str(index) + ".metastable.nc",
                                                      top=str(ARGS.topfile), image_molecules=bool(ARGS.periodic))
            index = index + 1
        if ARGS.ngl:
            VisualizeMetastable(my_samples, cmap, selection='backbone')


def Plot(FEATURE_TRAJ, FEATURES, CLUSTER, ITS, MSM, ARGS):
    # Concatenate data
    if ARGS.calc_tica_for_projection:
        feat_traj = FEATURE_TRAJ
    else:
        feat_traj = np.concatenate(FEATURE_TRAJ)

    # pyemma.plots.plot_feature_histograms(feat_traj, feature_labels=FEATURES);
    # plt.savefig(str(ARGS.output_prefix) + ".FeatureHistogram.png",dpi=300)

    # Plot free energy alone
    fig, axis = plt.subplots(1, 1)
    if ARGS.v_lim:
        min_value = int(ARGS.v_lim[0])
        max_value = int(ARGS.v_lim[1])
        pyemma.plots.plot_free_energy(*feat_traj.T[0:2], ax=axis, cmap=ARGS.cmap,
                                      vmin=min_value, vmax=max_value, legacy=False)
    else:
        pyemma.plots.plot_free_energy(*feat_traj.T[0:2], ax=axis, cmap=ARGS.cmap, legacy=False)
    axis.set_xlabel(ARGS.x_label)
    axis.set_ylabel(ARGS.y_label)
    if ARGS.x_lim:
        axis.set_xlim(float(ARGS.x_lim[0]), float(ARGS.x_lim[1]))
    if ARGS.y_lim:
        axis.set_ylim(float(ARGS.y_lim[0]), float(ARGS.y_lim[1]))
    fig.tight_layout()
    plt.savefig(str(ARGS.output_prefix) + ".FreeEnergy2D.png", dpi=300)
    plt.clf()

    # Plot corrected free energy alone
    fig, axis = plt.subplots(1, 1)
    if ARGS.v_lim:
        min_value = int(ARGS.v_lim[0])
        max_value = int(ARGS.v_lim[1])
        pyemma.plots.plot_free_energy(*feat_traj.T[0:2], weights=np.concatenate(MSM.trajectory_weights()), ax=axis,
                                      cmap=ARGS.cmap, vmin=min_value, vmax=max_value,
                                      legacy=False)  # bug if dim_red_stride and stride != 1
    else:
        pyemma.plots.plot_free_energy(*feat_traj.T[0:2], weights=np.concatenate(MSM.trajectory_weights()), ax=axis,
                              cmap=ARGS.cmap, legacy=False)  # bug if dim_red_stride and stride != 1
    axis.set_xlabel(ARGS.x_label)
    axis.set_ylabel(ARGS.y_label)
    if ARGS.x_lim:
        axis.set_xlim(float(ARGS.x_lim[0]), float(ARGS.x_lim[1]))
    if ARGS.y_lim:
        axis.set_ylim(float(ARGS.y_lim[0]), float(ARGS.y_lim[1]))
    fig.tight_layout()
    plt.savefig(str(ARGS.output_prefix) + ".FreeEnergy2D.Corrected.png", dpi=300)
    plt.clf()

    # Plot state diagram with kinetics
    if ARGS.n_metastable_states:
        # Count states
        nstates = int(ARGS.n_metastable_states)
        if not ARGS.hidden:
            MSM.pcca(nstates)
        print("Generating state map for " + str(nstates) + " markov model metastable states...")

        # Determine metastable states
        dtrajs_concatenated = np.concatenate(CLUSTER.dtrajs)
        metastable_traj = MSM.metastable_assignments[dtrajs_concatenated]

        # Plot state map
        fig, axis = plt.subplots(1, 1)
        _, _, misc = pyemma.plots.plot_state_map(*feat_traj.T[0:2], metastable_traj, ax=axis, legacy=False)
        misc['cbar'].set_ticklabels(range(0, nstates))

        # Add labels
        axis.set_xlabel(ARGS.x_label)
        axis.set_ylabel(ARGS.y_label)
        if ARGS.x_lim:
            axis.set_xlim(float(ARGS.x_lim[0]), float(ARGS.x_lim[1]))
        if ARGS.y_lim:
            axis.set_ylim(float(ARGS.y_lim[0]), float(ARGS.y_lim[1]))
        # Save and end
        fig.tight_layout()
        plt.savefig(str(ARGS.output_prefix) + ".StatesMap.png", dpi=300)
        plt.clf()

        # Make combined state map and corrected free energy plot
        fig, axis = plt.subplots(1, 2, sharex=True, sharey=True)
        pyemma.plots.plot_free_energy(*feat_traj.T[0:2], weights=np.concatenate(MSM.trajectory_weights()), ax=axis[0],
                                      cmap=ARGS.cmap, legacy=False)
        _, _, misc = pyemma.plots.plot_state_map(*feat_traj.T[0:2], metastable_traj, ax=axis[1], legacy=False)
        misc['cbar'].set_ticklabels(range(0, nstates))

        axis[0].set_xlabel(ARGS.x_label)
        axis[0].set_ylabel(ARGS.y_label)
        axis[1].set_xlabel(ARGS.x_label)
        axis[1].set_ylabel(ARGS.y_label)
        if ARGS.x_lim:
            axis[0].set_xlim(float(ARGS.x_lim[0]), float(ARGS.x_lim[1]))
            axis[1].set_xlim(float(ARGS.x_lim[0]), float(ARGS.x_lim[1]))
        if ARGS.y_lim:
            axis[0].set_ylim(float(ARGS.y_lim[0]), float(ARGS.y_lim[1]))
            axis[1].set_ylim(float(ARGS.y_lim[0]), float(ARGS.y_lim[1]))
        fig.tight_layout()
        plt.savefig(str(ARGS.output_prefix) + ".States.FreeEnergy2D.Corrected.png", dpi=300)
        plt.clf()

        highest_membership = MSM.metastable_distributions.argmax(1)
        coarse_state_centers = CLUSTER.clustercenters[MSM.active_set[int(nstates - 1)]]
        # coarse_state_centers = CLUSTER.clustercenters[MSM.active_set[highest_membership]]

        # Find pairwise
        # mfpt = np.zeros((nstates, nstates))
        # for i in range(nstates):
        #     for j in range(nstates):
        #         mfpt[i, j] = MSM.mfpt(
        #             MSM.metastable_sets[i],
        #             MSM.metastable_sets[j])
        #         print("MFPT states "+str(i)+" "+str(j)+" = "+str(mfpt[i, j]))
        #
        #
        # # Account for zeros
        # inverse_mfpt = np.zeros_like(mfpt)
        # nz = mfpt.nonzero()
        # inverse_mfpt[nz] = 1.0 / mfpt[nz]

        # Plot kinetic network
        # pyemma.plots.plot_network(
        #     inverse_mfpt,
        #     ax=axis,
        #     figpadding=0,
        #     pos=coarse_state_centers,
        #     arrow_label_format='%.1f ps',
        #     arrow_labels=mfpt,
        #     arrow_scale=1.0,
        #     state_scale=1.0,
        #     arrow_curvature=3.0,
        #     state_labels=range(0, nstates),
        #     show_frame=True,
        #     size=8)
        #
        # # Output kinetics figure
        # fig.tight_layout()
        # plt.savefig(str(ARGS.output_prefix) + ".FreeEnergy2D.States.png",dpi=300)
        # plt.clf()
        #
        # Create a committor / flux plot to show TPT
        if ARGS.hidden:
            A = [0]
            B = [1]
        else:
            A = MSM.metastable_sets[0]
            B = MSM.metastable_sets[1]
        flux = pyemma.msm.tpt(MSM, A, B)
        cg, cgflux = flux.coarse_grain(MSM.metastable_sets)

        if ARGS.hidden:
            fig, ax = plt.subplots(1, 1)
            pyemma.plots.plot_contour(
                *feat_traj.T,
                metastable_traj,
                #cmap='brg',
                #ax=ax,
                mask=True,
                cbar_label=r'committor A $\to$ B',
                alpha=0.8,
                zorder=-1);
        else:
            fig, ax = plt.subplots(1, 1)
            pyemma.plots.plot_contour(
                *feat_traj.T,
                flux.committor[dtrajs_concatenated],
                #cmap='brg',
                #ax=ax,
                mask=True,
                cbar_label=r'committor A $\to$ B',
                alpha=0.8,
                zorder=-1);

        pyemma.plots.plot_flux(
            cgflux,
            coarse_state_centers,
            cgflux.stationary_distribution,
            state_labels=['A','B'],
            ax=ax,
            show_committor=False,
            figpadding=0,
            show_frame=True,
            arrow_label_format='%2.e / ps');

        # Add labels
        ax.set_xlabel(ARGS.x_label)
        ax.set_ylabel(ARGS.y_label)
        if ARGS.x_lim:
            ax.set_xlim(float(ARGS.x_lim[0]), float(ARGS.x_lim[1]))
        if ARGS.y_lim:
            ax.set_ylim(float(ARGS.y_lim[0]), float(ARGS.y_lim[1]))

        fig.tight_layout()
        plt.savefig(str(ARGS.output_prefix) + ".FreeEnergy2D.TPT.png",dpi=300)
        plt.clf()

        # paths, path_fluxes = cgflux.pathways(fraction=0.99)
        # print('percentage       \tpath')
        # print('-------------------------------------')
        # for i in range(len(paths)):
        #     print(np.round(path_fluxes[i] / np.sum(path_fluxes), 3),' \t', paths[i] + 1)

    # Create a combined plot
    if len(ARGS.lag) > 1:
        figc, axes = plt.subplots(1, 3, figsize=(12, 7))
        pyemma.plots.plot_feature_histograms(feat_traj, feature_labels=[str(ARGS.x_label), str(ARGS.y_label)],
                                             ax=axes[0])
        pyemma.plots.plot_density(*feat_traj.T, ax=axes[1], cbar=False, alpha=0.1)
        axes[1].scatter(*CLUSTER.clustercenters.T, s=2, c='C1')
        axes[1].set_xlabel(ARGS.x_label)
        axes[1].set_ylabel(ARGS.y_label)
        if ARGS.x_lim:
            axes[1].set_xlim(float(ARGS.x_lim[0]), float(ARGS.x_lim[1]))
        if ARGS.y_lim:
            axes[1].set_ylim(float(ARGS.y_lim[0]), float(ARGS.y_lim[1]))
        axes[1].set_aspect('equal')
        pyemma.plots.plot_implied_timescales(ITS, ylog=True, ax=axes[2])
        figc.tight_layout()
        plt.savefig(str(ARGS.output_prefix) + ".Combined.png", dpi=300)


def VisualizeMetastable(samples, cmap, selection='backbone'):
    """ visualize metastable states
    Parameters
    ----------
    samples: list of mdtraj.Trajectory objects
        each element contains all samples for one metastable state.
    cmap: matplotlib.colors.ListedColormap
        color map used to visualize metastable states before.
    selection: str
        which part of the molecule to selection for visualization. For details have a look here:
        http://mdtraj.org/latest/examples/atom-selection.html#Atom-Selection-Language
    """
    import nglview
    from matplotlib.colors import to_hex

    widget = nglview.NGLWidget()
    widget.clear_representations()
    ref = samples[0]
    for i, s in enumerate(samples):
        s = s.superpose(ref)
        s = s.atom_slice(s.top.select(selection))
        comp = widget.add_trajectory(s)
        comp.add_ball_and_stick()

    # this has to be done in a separate loop for whatever reason...
    x = np.linspace(0, 1, num=len(samples))
    for i, x_ in enumerate(x):
        c = to_hex(cmap(x_))
        widget.update_ball_and_stick(color=c, component=i, repr_index=i)
        widget.remove_cartoon(component=i)
    return widget


# Argument parser
def cmdlineparse():
    # Construct argument parser
    parser = ArgumentParser(description="Command line arguments", add_help=True)

    # Make multiple argument groups for organization purposes
    global_group = parser.add_argument_group("Global", "General arguments")
    feature_group = parser.add_argument_group("Features",
                                              "Indicate what types of features to generate for the trajectory")
    ref_feature_group = parser.add_argument_group("Reference Features",
                                                  "Indicate features to be processed with TICA onto which input features will be projected")
    dimensions_group = parser.add_argument_group("Dimensionality Reduction and Clustering",
                                                 "Arguments related to dimensionality reduction and clustering methods")
    plot_group = parser.add_argument_group("Plotting", "Specify plotting details (e.g. axis label names)")

    # Global group
    global_group.add_argument("-topfile", dest="topfile", required=True, help="Topology file",
                              metavar="<topology file>")
    global_group.add_argument("-trajfiles", dest="trajfiles", nargs="+", required=False, help="Trajectory files",
                              metavar="<trajectory files>")
    global_group.add_argument("-output_prefix", dest="output_prefix", required=False, default="PyemmaOutput",
                              help="Output prefix for output files", metavar="<output files>")
    global_group.add_argument("-stride", dest="stride", required=False, default=1,
                              help="Load every <stride>th trajectory frame", metavar="<stride>")
    global_group.add_argument("-dt", dest="dt", required=False, default='1 step', help="Frame step quantity and unit",
                              metavar="<dt>")
    global_group.add_argument("-lag", dest="lag", required=False, default='1', nargs="+", help="Number of lag steps",
                              metavar="<lag>")
    global_group.add_argument("-periodic", dest="periodic", required=False, default=False, action="store_true",
                              help="Indicate the presence of periodic boundary conditions")
    global_group.add_argument("-checkpoint_save", dest="checkpoint_save", required=False, nargs=int(2),
                              help="Specify a file name and model name to save the clustered features externally for faster loading or mapping",
                              metavar="<file_name model_name>")
    global_group.add_argument("-checkpoint_load", dest="checkpoint_load", required=False, nargs=int(2),
                              help="Specify a file name and model name base (i.e. without suffix) to load all data for that model",
                              metavar="<file_name model_name>")
    global_group.add_argument("-load_tica_for_projection", dest="load_tica_for_projection", required=False,
                              help="Specify a file name and model name base (i.e. without suffix) to load a tica featurized trajectory to project input trajectory onto",
                              metavar="<file_name model_name>")
    global_group.add_argument("-calc_tica_for_projection", dest="calc_tica_for_projection", required=False, nargs="+",
                              help="specify a secondary topfile and trajectory, compute tica, and load main input data onto this tica trajectory",
                              metavar="<topfile trajectory>")
    global_group.add_argument("-hidden", dest="hidden", required=False, default=False, action="store_true",
                              help="Generate HMM instead of regular MSM")

    # Feature group
    feature_group.add_argument("-residue_pair_sc_files", dest="residue_pair_sc_files", required=False, nargs="+",
                               help="Compute the contact distance between nearest sidechain heavy atoms in these pairs of residues",
                               metavar="<residue_pair_sc_files>")
    feature_group.add_argument("-residue_pair_closest_files", dest="residue_pair_closest_files", required=False,
                               nargs="+",
                               help="Compute the contact distance between the closest heavy atoms in these pairs of residues",
                               metavar="<residue_pair_closest_files>")
    feature_group.add_argument("-dihedral_atoms_files", dest="dihedral_atoms_files", required=False, nargs="+",
                               help="Compute the dihedral angles formed by these sets of four atoms",
                               metavar="<dihedral_atoms_files>")
    feature_group.add_argument("-atom_pairs_files", dest="atom_pairs_files", required=False, nargs="+",
                               help="Compute the distances formed by each pair of atoms", metavar="<atom_pairs_files>")
    feature_group.add_argument("-atom_contacts_files", dest="atom_contacts_files", required=False, nargs="+",
                               help="Computes whether contacts are formed by each pair of atoms",
                               metavar="<atom_contacts_files>")
    feature_group.add_argument("-atom_groups_files", dest="atom_groups_files", required=False, nargs="+",
                               help="Compute distances between groups specified here based on group_pairs",
                               metavar="<atom_groups_files>")
    feature_group.add_argument("-atom_iterables_a_files", dest="atom_iterables_a_files", required=False, nargs="+",
                               help="Compute distances between atoms specified here or those specified in atom_iterables_b",
                               metavar="<atom_iterables_a>")
    feature_group.add_argument("-atom_iterables_b_files", dest="atom_iterables_b_files", required=False, nargs="+",
                               help="Compute distances between atoms specified in atom_iterables_a",
                               metavar="<atom_iterables_b>")
    feature_group.add_argument("-group_pairs_file", dest="group_pairs_file", required=False,
                               help="Which groups to compute distances between", metavar="<group_pairs_file>")
    feature_group.add_argument("-backbone_torsions", dest="backbone_torsions", required=False, default="",
                               help="Add backbone torsions to the feature list")
    feature_group.add_argument("-sidechain_torsions", dest="sidechain_torsions", required=False, default="",
                               help="Add sidechain torsions to the feature list")

    # Reference feature group
    ref_feature_group.add_argument("-ref_residue_pair_sc_files", dest="ref_residue_pair_sc_files", required=False,
                                   nargs="+",
                                   help="Compute the contact distance between nearest sidechain heavy atoms in these pairs of residues",
                                   metavar="<residue_pair_sc_files>")
    ref_feature_group.add_argument("-ref_residue_pair_closest_files", dest="ref_residue_pair_closest_files",
                                   required=False, nargs="+",
                                   help="Compute the contact distance between the closest heavy atoms in these pairs of residues",
                                   metavar="<residue_pair_closest_files>")
    ref_feature_group.add_argument("-ref_dihedral_atoms_files", dest="ref_dihedral_atoms_files", required=False,
                                   nargs="+",
                                   help="Compute the dihedral angles formed by these sets of four atoms",
                                   metavar="<dihedral_atoms_files>")
    ref_feature_group.add_argument("-ref_atom_pairs_files", dest="ref_atom_pairs_files", required=False, nargs="+",
                                   help="Compute the distances formed by each pair of atoms",
                                   metavar="<atom_pairs_files>")
    ref_feature_group.add_argument("-ref_atom_groups_files", dest="ref_atom_groups_files", required=False, nargs="+",
                                   help="Compute distances between groups specified here based on group_pairs",
                                   metavar="<atom_groups_files>")
    ref_feature_group.add_argument("-ref_group_pairs_file", dest="ref_group_pairs_file", required=False,
                                   help="Which groups to compute distances between", metavar="<group_pairs_file>")
    ref_feature_group.add_argument("-ref_backbone_torsions", dest="ref_backbone_torsions", required=False,
                                   action="store_true", default=False,
                                   help="Add all backbone torsions to the feature list")
    ref_feature_group.add_argument("-ref_sidechain_torsions", dest="ref_sidechain_torsions", required=False,
                                   action="store_true", default=False,
                                   help="Add all sidechain torsions to the feature list")

    # Dimensions group
    dimensions_group.add_argument("-tica", dest="tica", required=False,
                                  help="Perform TICA with the requested kinetic variance cutoff",
                                  metavar="<tica_cutoff>")
    dimensions_group.add_argument("-nits", dest="nits", required=False, default=5,
                                  help="Number of implied timescales",
                                  metavar="nits")
    dimensions_group.add_argument("-ck_test", dest="ck_test", required=False, default=False,
                                  help="Perform Chapman-Kolmogorow test with the specified number of metastable states")
    dimensions_group.add_argument("-pca", dest="pca", required=False,
                                  help="Perform PCA with the requested variance explained",
                                  metavar="<variance_explained>")
    dimensions_group.add_argument("-kmeans", dest="kmeans", required=False, default=False,
                                  help="Perform kmeans clustering of dimensionality-reduced features (good if your data has several well-sampled states)",
                                  action="store_true")
    dimensions_group.add_argument("-k_clusters", dest="k_clusters", required=False, default=100,
                                  help="Number of clusters for k-means clustering",
                                  metavar="<k_klusters>")
    dimensions_group.add_argument("-k_iter", dest="k_iter", required=False, default=10,
                                  help="Number of iterations for k-means clustering",
                                  metavar="<k_iter>")
    dimensions_group.add_argument("-regspace", dest="regspace", required=False, default=False,
                                  help="Perform regular space clustering of dimensionality-reduced features (good if you want to represent sparsely sampled states, too)",
                                  action="store_true")
    dimensions_group.add_argument("-regspace_dmin", dest="regspace_dmin", required=False, default=0.50,
                                  help="Dmin spacing for regspace clustering",
                                  metavar="<regspace_dmin>")
    dimensions_group.add_argument("-dim_red_stride", dest="dim_red_stride", required=False, default=1,
                                  help="Load every <stride>th trajectory frame for dimensionality reduction tasks; distinct from trajectory stride",
                                  metavar="<stride>")
    dimensions_group.add_argument("-n_metastable_states", dest="n_metastable_states", required=False,
                                  help="Number of metastable states for network analysis",
                                  metavar="<n_metastable_states>")
    dimensions_group.add_argument("-n_samples_each_state", dest="n_samples_each_state", required=False, default=1000,
                                  help="Number of structures to sample from each metastable state",
                                  metavar="<n_samples_each_state>")

    # Plots group
    plot_group.add_argument("-x_label", dest="x_label", required=False, default="X",
                            help="X axis label for free energy plot", metavar="<x label>")
    plot_group.add_argument("-y_label", dest="y_label", required=False, default="Y",
                            help="Y axis label for free energy plot", metavar="<y label>")
    plot_group.add_argument("-x_lim", dest="x_lim", required=False,
                            help="X axis min and max values for free energy plots", nargs="+", metavar="<xmin xmax>")
    plot_group.add_argument("-y_lim", dest="y_lim", required=False,
                            help="Y axis min and max values for free energy plots", nargs="+", metavar="<ymin ymax>")
    plot_group.add_argument("-v_lim", dest="v_lim", required=False,
                            help="Colorbar axis min and max values for free energy plots", nargs="+", metavar="<ymin ymax>")
    plot_group.add_argument("-ngl", dest="ngl", required=False, default=False, action="store_true",
                            help="Visuable metastable states in NGL viewer")
    plot_group.add_argument("-cmap", dest="cmap", required=False, default="nipy_spectral",
                            help="cmap string for free energy plots")
    plot_group.add_argument("-states_cmap", dest="states_cmap", required=False, default="Set3",
                            help="cmap string for state map plots")

    # Done
    args = parser.parse_args()
    return args


# Run from commandline
if __name__ == '__main__':
    main()
