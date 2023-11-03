    n_vert = len(load(surface))

    # get data from dlabel / medial wall files if provided
    labels, mask = None, np.zeros(n_vert, dtype=bool)
    if dlabel is not None:
        labels = check_image_file(dlabel)
        if len(labels) != n_vert:
            raise ValueError('Provided `dlabel` file does not contain same '
                             'number of vertices as provided `surface`')
    if medial is not None:
        mask = np.asarray(check_image_file(medial), dtype=bool)
        if len(mask) != n_vert:
            raise ValueError('Provided `medial` file does not contain same '
                             'number of vertices as provided `surface`')

    # define which function we'll be using to calculate the distances
    if euclid:
        func = _get_euclid_distance
        graph = check_surface(surface)  # vertex coordinates
    else:
        if use_wb:
            func = _get_workbench_distance
            graph = surface
        else:
            func = _get_graph_distance
            vert, faces = [darray.data for darray in nib.load(surface).darrays]
            graph = make_surf_graph(vert, faces, mask=mask)

    # if we want the vertex-vertex distance matrix we'll stream it to disk to
    # save on memory, a la `_geodesic()` or `_euclid()`
    # NOTE: streaming to disk takes a lot more _time_ than storing in memory
    if labels is None:
        with open(outfile, 'w') as dest:
            for n in range(n_vert):
                if verbose and n % 1000 == 0:
                    print('Running vertex {} of {}'.format(n, n_vert))
                np.savetxt(dest, func(n, graph))
    # we can store the temporary n_vert x label matrix in memory; running this
    # is much faster than trying to read through the giant vertex-vertex
    # distance matrix file
    else:
        # depends on size of parcellation, but assuming even a liberal 1000
        # parcel atlas this will be ~250 MB in-memory for the default fslr32k
        # resolution
        unique_parcels = np.unique(labels)
        par, func = Parallel(n_jobs=n_jobs), delayed(func)
        dist = np.row_stack(par(func(n, graph, labels) for n in range(n_vert)))
        # average rows (vertices) into parcels; columns are already parcels
        dist = np.row_stack([
            dist[labels == lab].mean(axis=0) for lab in unique_parcels])
        dist[np.diag_indices_from(dist)] = 0
        # NOTE: if `medial` is supplied and any of the parcel labels correspond
        # to the medial wall then those parcel-parcel distances will be `inf`!

        # remove unassigned parcel
        if unassigned_value in unique_parcels:
            idx = list(unique_parcels).index(unassigned_value)
            dist = np.delete(dist, idx, axis=0)
            dist = np.delete(dist, idx, axis=1)

        np.savetxt(outfile, dist)

    return outfile