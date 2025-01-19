#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <math.h>
#include <cstdlib>
#include "Graph.h"
#include "Problem.h"

using namespace std;
const vector<string> DSJ_FILES = {
    "DSJC1000.1.col",
    "DSJC1000.5.col",
    "DSJC1000.9.col",
    "DSJC125.1.col",
    "DSJC125.5.col",
    "DSJC125.9.col",
    "DSJC250.1.col",
    "DSJC250.5.col",
    "DSJC250.9.col",
    "DSJC500.1.col",
    "DSJC500.5.col",
    "DSJC500.9.col",
    "DSJR500.1.col",
    "DSJR500.1c.col",
    "DSJR500.5.col",
};
const vector<string> REG_FILES = {
    "fpsol2.i.1.col",
    "fpsol2.i.2.col",
    "fpsol2.i.3.col",
    "inithx.i.1.col",
    "inithx.i.2.col",
    "inithx.i.3.col",
    "mulsol.i.1.col",
    "mulsol.i.2.col",
    "mulsol.i.3.col",
    "mulsol.i.4.col",
    "mulsol.i.5.col",
    "zeroin.i.1.col",
    "zeroin.i.2.col",
    "zeroin.i.3.col",
};
const vector<string> LAT_FILES = {
    "latin_square_10.col",
};
const vector<string> LEI_FILES = {
    "le450_15a.col",
    "le450_15b.col",
    "le450_15b.col.txt",
    "le450_15c.col",
    "le450_15d.col",
    "le450_25a.col",
    "le450_25b.col",
    "le450_25c.col",
    "le450_25d.col",
    "le450_5a.col",
    "le450_5b.col",
    "le450_5c.col",
    "le450_5d.col",
};
const vector<string> SCH_FILES = {
    "school1.col",
    "school1_nsh.col",
};
const vector<string> SGB_FILES = {
    "anna.col",
    "david.col",
    "games120.col",
    "homer.col",
    "huck.col",
    "jean.col",
    "miles1000.col",
    "miles1500.col",
    "miles250.col",
    "miles500.col",
    "miles750.col",
    "queen10_10.col",
    "queen11_11.col",
    "queen12_12.col",
    "queen13_13.col",
    "queen14_14.col",
    "queen15_15.col",
    "queen16_16.col",
    "queen5_5.col",
    "queen6_6.col",
    "queen7_7.col",
    "queen8_12.col",
    "queen8_8.col",
    "queen9_9.col",
};
const vector<string> MYC_FILES = {
    "myciel3.col",
    "myciel4.col",
    "myciel5.col",
    "myciel6.col",
    "myciel7.col",
};
const vector<string> MIZ_FILES = {
    "mug100_1.col",
    "mug100_25.col",
    "mug88_1.col",
    "mug88_25.col",
};
const vector<string> HOS_FILES = {};
const vector<string> CAR_FILES = {
    "1-FullIns_3.col",
    "1-FullIns_4.col",
    "1-FullIns_5.col",
    "1-Insertions_4.col",
    "1-Insertions_5.col",
    "1-Insertions_6.col",
    "2-FullIns_3.col",
    "2-FullIns_4.col",
    "2-FullIns_5.col",
    "2-Insertions_3.col",
    "2-Insertions_4.col",
    "2-Insertions_5.col",
    "3-FullIns_3.col",
    "3-FullIns_4.col",
    "3-FullIns_5.col",
    "3-Insertions_3.col",
    "3-Insertions_4.col",
    "3-Insertions_5.col",
    "4-FullIns_3.col",
    "4-FullIns_4.col",
    "4-FullIns_5.col",
    "4-Insertions_3.col",
    "4-Insertions_4.col",
    "5-FullIns_3.col",
    "5-FullIns_4.col",
};

struct DataGraph
{
    string nameGraph;
    HeuristicaConstrutivaResponse data;
};

void processaInstancia(vector<string> namesFiles, vector<Graph> &instancesGraphs)
{
    Problem problem;
    for (int i = 0; i < namesFiles.size(); i++)
    {
        Graph graph;
        graph.read_from_file(namesFiles[i]);
        instancesGraphs.push_back(graph);
        HeuristicaConstrutivaResponse response = problem.heuristicaConstrutiva(graph);
        cout << "-- Grafo [" << namesFiles[i] << "]com nVertices: " << graph.get_num_vertices();
        cout << " - Cores: " << response.num_cores_encontradas
             << " - Clique: " << response.size_clique
             << " - Delta: " << response.delta_clique
             << " - Arestas: " << response.num_arestas
             << " - DArestas: " << response.densidade_arestas << endl;
    }
};

int main()
{

    ifstream file("grafo_aleatorio.txt");
    if (!file.is_open())
    {
        cout << "Error opening file" << endl;
        return 1;
    }

    string line;
    // vector<Graph> randomGraphs;
    // int actual_graph = -1;
    // while (getline(file, line))
    // {
    //     stringstream linestream;
    //     linestream << line;
    //     char info[10];
    //     linestream.getline(info, 4, ':');
    //     if (!line.empty() && atoi(info) != actual_graph)
    //     {
    //         if (actual_graph != -1)
    //         {
    //             randomGraphs[actual_graph].generate_adj_matrix();
    //             cout << "-- Grafo [" << i << "]com nVertices: " << randomGraphs[i].get_num_vertices();
    //             HeuristicaConstrutivaResponse response = problem.heuristicaConstrutiva(randomGraphs[i]);
    //             cout << " - Cores: " << response.num_cores_encontradas
    //                  << " - Clique: " << response.size_clique
    //                  << " - Delta: " << response.delta_clique
    //                  << " - Arestas: " << response.num_arestas
    //                  << " - DArestas: " << response.densidade_arestas << endl;
    //         }
    //         actual_graph++;
    //         Graph graph;
    //         randomGraphs.push_back(graph);
    //     }

    //     randomGraphs[actual_graph].get_from_line(line);
    // }
    // file.close();
    // cout << "Graphs: " << randomGraphs.size() << endl;

    // mostraGrafo(graphs[0]);

    Problem problem;

    // int initial_k = 5;

    // for (int i = 0; i < randomGraphs.size(); i++)
    // {
    //     if (randomGraphs[i].get_num_vertices() <= 300)
    //     {
    //         continue;
    //     }
    //     cout << "-- Grafo [" << i << "]com nVertices: " << randomGraphs[i].get_num_vertices() << flush;
    //     int k = initial_k;
    //     cout << " tentativas de k: " << flush;
    //     while (problem.temKCliqueBK(randomGraphs[i], k))
    //     {
    //         cout << k << "," << flush;
    //         k++;
    //     }
    //     cout << endl;
    // }

    // int media_100_vertices = 0;
    // float media_delta_100_vertices = 0;
    // int count_100_vertices = 0;
    // int media_300_vertices = 0;
    // float media_delta_300_vertices = 0;
    // int count_300_vertices = 0;
    // int media_500_vertices = 0;
    // float media_delta_500_vertices = 0;
    // int count_500_vertices = 0;
    // int media_1000_vertices = 0;
    // float media_delta_1000_vertices = 0;
    // int count_1000_vertices = 0;

    // for (int i = 0; i < randomGraphs.size(); i++)
    // {
    //     cout << "-- Grafo [" << i << "]com nVertices: " << randomGraphs[i].get_num_vertices();
    //     HeuristicaConstrutivaResponse response = problem.heuristicaConstrutiva(randomGraphs[i]);
    //     cout << " - Cores: " << response.num_cores_encontradas
    //          << " - Clique: " << response.size_clique
    //          << " - Delta: " << response.delta_clique
    //          << " - Arestas: " << response.num_arestas
    //          << " - DArestas: " << response.densidade_arestas << endl;

    //     switch (randomGraphs[i].get_num_vertices())
    //     {
    //     case 100:
    //         media_100_vertices += response.num_cores_encontradas;
    //         media_delta_100_vertices += response.delta_clique;
    //         count_100_vertices++;
    //         break;
    //     case 300:
    //         media_300_vertices += response.num_cores_encontradas;
    //         media_delta_300_vertices += response.delta_clique;
    //         count_300_vertices++;
    //         break;
    //     case 500:
    //         media_500_vertices += response.num_cores_encontradas;
    //         media_delta_500_vertices += response.delta_clique;
    //         count_500_vertices++;
    //         break;
    //     case 1000:
    //         media_1000_vertices += response.num_cores_encontradas;
    //         media_delta_1000_vertices += response.delta_clique;
    //         count_1000_vertices++;
    //         break;
    //     default:
    //         break;
    //     }
    // }

    // media_100_vertices /= count_100_vertices;
    // media_delta_100_vertices /= count_100_vertices;
    // media_300_vertices /= count_300_vertices;
    // media_delta_300_vertices /= count_300_vertices;
    // media_500_vertices /= count_500_vertices;
    // media_delta_500_vertices /= count_500_vertices;
    // media_1000_vertices /= count_1000_vertices;
    // media_delta_1000_vertices /= count_1000_vertices;

    // cout << "Media 100 vertices: " << media_100_vertices;
    // cout << " Media delta: " << media_100_vertices << endl;
    // cout << "Media 300 vertices: " << media_300_vertices;
    // cout << " Media delta: " << media_300_vertices << endl;
    // cout << "Media 500 vertices: " << media_500_vertices;
    // cout << " Media delta: " << media_500_vertices << endl;
    // cout << "Media 1000 vertices: " << media_1000_vertices;
    // cout << " Media delta: " << media_1000_vertices << endl;

    // for (int i = 16; i < 20; i++)
    // {
    //     cout << "-----------------------------" << endl;
    //     cout << "1 rodada: " << i << endl;
    //     cout << "Grafo com nVertices: " << graphs[0].get_num_vertices() << endl;
    //     Problem problem;
    //     int n_cores = i;
    //     int max_iter = 1000;

    //     // if (problem.tabucol(graphs[0], max_iter, n_cores))
    //     // {
    //     //     break;
    //     // }

    //     if (problem.tabucol2(graphs[0], max_iter, n_cores))
    //     {
    //         break;
    //     }
    // }

    // for (int i = 35; i < 50; i++)
    // {
    //     cout << "-----------------------------" << endl;
    //     cout << "1 rodada: " << i << endl;
    //     cout << "Grafo com nVertices: " << graphs[20].get_num_vertices() << endl;
    //     Problem problem;
    //     int n_cores = i;
    //     int max_iter = 20000;

    //     // if (problem.tabucol(graphs[0], max_iter, n_cores))
    //     // {
    //     //     break;
    //     // }

    //     if (problem.tabucol2(graphs[20], max_iter, n_cores))
    //     {
    //         break;
    //     }
    // }

    // for (int i = 50; i < 75; i++)
    // {
    //     cout << "-----------------------------" << endl;
    //     cout << "1 rodada: " << i << endl;
    //     cout << "Grafo com nVertices: " << graphs[31].get_num_vertices() << endl;
    //     Problem problem;
    //     int n_cores = i;
    //     int max_iter = 1000000;

    //     // if (problem.tabucol(graphs[0], max_iter, n_cores))
    //     // {
    //     //     break;
    //     // }

    //     if (problem.tabucol2(graphs[31], max_iter, n_cores))
    //     {
    //         break;
    //     }
    // }

    // for (int i = 85; i < 100; i++)
    // {
    //     cout << "-----------------------------" << endl;
    //     cout << "1 rodada: " << i << endl;
    //     cout << "Grafo com nVertices: " << graphs[36].get_num_vertices() << endl;
    //     Problem problem;
    //     int n_cores = i;
    //     int max_iter = 1000000;

    //     // if (problem.tabucol(graphs[0], max_iter, n_cores))
    //     // {
    //     //     break;
    //     // }

    //     if (problem.tabucol2(graphs[36], max_iter, n_cores))
    //     {
    //         break;
    //     }
    // }

    // for (int i = 0; i < graphs.size(); i++)
    // {
    //     cout << "Graph " << i << " vertices: " << graphs[i].get_num_vertices() << endl;
    // }

    // vector<string> namesFiles = {
    //     "le450_15b.col.txt",
    // };

    int actual_graph = -1;
    Graph *graph;
    while (getline(file, line))
    {
        stringstream linestream;
        linestream << line;
        char info[10];
        linestream.getline(info, 4, ':');
        if (!line.empty() && atoi(info) != actual_graph)
        {
            if (actual_graph != -1)
            {
                graph->generate_adj_matrix();
                cout << "-- Grafo [" << actual_graph << "]com nVertices: " << graph->get_num_vertices();
                HeuristicaConstrutivaResponse response = problem.heuristicaConstrutiva(*graph);
                cout << " - Cores: " << response.num_cores_encontradas
                     << " - Clique: " << response.size_clique
                     << " - Delta: " << response.delta_clique
                     << " - Arestas: " << response.num_arestas
                     << " - DArestas: " << response.densidade_arestas << endl;
            }
            actual_graph++;
            graph = new Graph();
        }

        graph->get_from_line(line);
    }
    file.close();
    vector<Graph> instancesGraphs;

    processaInstancia(DSJ_FILES, instancesGraphs);
    processaInstancia(REG_FILES, instancesGraphs);
    processaInstancia(LAT_FILES, instancesGraphs);
    processaInstancia(LEI_FILES, instancesGraphs);
    processaInstancia(SCH_FILES, instancesGraphs);
    processaInstancia(SGB_FILES, instancesGraphs);
    processaInstancia(MYC_FILES, instancesGraphs);
    processaInstancia(MIZ_FILES, instancesGraphs);
    processaInstancia(HOS_FILES, instancesGraphs);
    processaInstancia(CAR_FILES, instancesGraphs);

    // for (int i = 0; i < namesFiles.size(); i++)
    // {
    //     Graph graph;
    //     graph.read_from_file(namesFiles[i]);
    //     instancesGraphs.push_back(graph);
    //     HeuristicaConstrutivaResponse response = problem.heuristicaConstrutiva(graph);
    //     cout << " - Cores: " << response.num_cores_encontradas
    //          << " - Grau maximo: " << response.num_grau_maximo
    //          << " - Delta: " << response.delta
    //          << " - Arestas: " << response.num_arestas
    //          << " - DArestas: " << response.densidade_arestas << endl;
    // }

    return 0;
}
