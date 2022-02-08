int mpi_vertex_dist(graph_t *graph, int start_vertex, int *result)
{
    int num_vertices = graph->num_vertices; 
    fill_n(result, num_vertices, MAX_DIST);
    
    auto start_time = Time::now();
    
    int depth = 0;
    result[start_vertex] = depth;

    int keep_going = true;

    int rank = MPI::COMM_WORLD.Get_rank();
    int size = MPI::COMM_WORLD.Get_size();
    int neighbors_vertices[num_vertices] = {0};
    int neighbors_length = 0;
    int tag0 = 0;
    int tag1 = 1;

    while (keep_going)
    {
        keep_going = false;

        for (int vertex = 0; vertex < num_vertices; vertex++) {

            if (vertex % size != rank) continue;

            if (result[vertex] == depth) {
                for (int n = graph->v_adj_begin[vertex]; 
                    n < graph->v_adj_begin[vertex] + graph->v_adj_length[vertex]; 
                    n++)
                {
                    int neighbor = graph->v_adj_list[n];

                    if (result[neighbor] > depth+1)
                    {
                        if (rank == 0) {
                            result[neighbor] = depth+1;
                            keep_going = true;
                        } else {
                            neighbors_vertices[neighbors_length] = neighbor;
                            neighbors_length++;
                        }
                    }
                }
            }
        }
        if (rank == 0) {
            for (int i=1;i<size;i++) {
                MPI_Recv(&neighbors_length,1,MPI_INT,i,tag0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                MPI_Recv(neighbors_vertices,neighbors_length,MPI_INT,i,tag1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                for (int j=0;j<neighbors_length;j++) {
                    keep_going = true;
                    int neighbor_index = neighbors_vertices[j];
                    result[neighbor_index] = depth+1;
                }
            }
        } else {
            MPI_Send(&neighbors_length,1,MPI_INT,0,tag0,MPI_COMM_WORLD);
            MPI_Send(neighbors_vertices,neighbors_length,MPI_INT,0,tag1,MPI_COMM_WORLD);
        }
        MPI_Bcast(&keep_going, 1, MPI_CXX_BOOL, 0, MPI_COMM_WORLD);
        MPI_Bcast(result, num_vertices, MPI_INT, 0, MPI_COMM_WORLD);
        neighbors_length = 0;
        depth++;
    }
    return std::chrono::duration_cast<us>(Time::now()-start_time).count();
}