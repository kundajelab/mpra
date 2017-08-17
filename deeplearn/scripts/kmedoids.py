import numpy as np


def _kmeanspp(dist, k):
    num_points = dist.shape[0]
    centers = np.empty(k, dtype=np.int)
    centers[0] = np.random.randint(num_points)

    for i in xrange(1, k):
        # get closest points
        closest_center = np.argmin(dist[:, centers[:i]], axis=1)
        dist_sq = \
            np.square(dist[np.arange(num_points), centers[closest_center]])
        next_center = np.random.choice(np.arange(num_points),
                                       p=dist_sq / np.sum(dist_sq))
        centers[i] = next_center

    return centers


def kmedoids(dist, n_clusters, n_init=5, max_iter=300, init='kmeans++',
             random_search_tries=3):
    best_cost = np.inf

    for _ in xrange(n_init):
        cluster_idx, medoids, cost = \
            _kmedoids(dist, n_clusters, max_iter, init,
                      random_search_tries)

        if cost < best_cost:
            best_cost = cost
            best_vals = (cluster_idx, medoids, cost)

    return best_vals


def _reassign_clusters(dist, medoids):
    num_points = dist.shape[0]

    cluster_idx = np.argmin(dist[:, medoids], axis=1)
    cost = dist[np.arange(num_points), medoids[cluster_idx]].sum()
    return cluster_idx, cost


def _kmedoids(dist, k, max_iter, init, num_tries=3,
              max_iters_without_improvement=1):
    # dist should be a symmetric matrix of distances
    assert dist.ndim == 2 and dist.shape[0] == dist.shape[1]

    num_points = dist.shape[0]

    if init == 'random':
        medoids = np.random.choice(num_points, k, replace=False)
    elif init == 'kmeans++':
        medoids = _kmeanspp(dist, k)
    else:
        raise ValueError('Invalid initalization.')

    min_cost = np.inf
    best_medoids = np.copy(medoids)
    iters_without_improvement = 0

    cluster_idx, cost = _reassign_clusters(dist, medoids)
    print('Initial cost: {}'.format(cost))

    for iter_ in xrange(max_iter):
        # recompute medoids
        for medoid_idx in xrange(k):
            points_in_cluster = np.where(cluster_idx == medoid_idx)[0]
            submatrix_index = np.ix_(points_in_cluster, points_in_cluster)
            dist_sums = dist[submatrix_index].sum(axis=1)
            new_medoid = np.argmin(dist_sums)
            medoids[medoid_idx] = points_in_cluster[new_medoid]

        cluster_idx, cost = _reassign_clusters(dist, medoids)

        # monte carlo search
        for medoid_idx in xrange(k):
            improvement = False
            original_medoid = medoids[medoid_idx]

            proposal_idx = np.random.randint(0, num_points, size=num_tries)
            
            new_best_idx = 0
            for proposal in proposal_idx:
                medoids[medoid_idx] = proposal
                new_cluster_idx, new_cost = _reassign_clusters(dist, medoids)
                if new_cost < cost:
                    if new_best_idx % 15 == 0:
                        print('Random search cost: {}'.format(new_cost))
                    new_best_idx += 1
                    cost = new_cost
                    best_cluster_idx = new_cluster_idx
                    best_reassignment = proposal
                    improvement = True

            if improvement:
                medoids[medoid_idx] = best_reassignment
                cluster_idx = best_cluster_idx
            else:
                medoids[medoid_idx] = original_medoid

        if cost < min_cost:
            min_cost = cost
            iters_without_improvement = 0
            best_medoids[:] = medoids
        else:
            iters_without_improvement += 1

        print('Iter {} - Cost: {}, Min cost: {}'.format(iter_, cost, min_cost))

        if iters_without_improvement >= max_iters_without_improvement:
            print('Terminating: no improvement in the last {} steps'.format(
                iters_without_improvement))
            break

    cluster_idx = np.argmin(dist[:, best_medoids], axis=1)
    return cluster_idx, best_medoids, min_cost
