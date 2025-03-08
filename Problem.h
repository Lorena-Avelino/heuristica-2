#ifndef PROBLEM_H
#define PROBLEM_H

#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <limits>
#include <unordered_map>
#include <chrono>
#include <set>
#include <algorithm>
#include <unordered_set>
#include "Graph.h"
#include <deque>
#include <queue>
#include <math.h>

using namespace std;
using namespace chrono;

struct HeuristicaConstrutivaResponse
{
    int num_cores_encontradas;
    int num_vertices;
    int num_grau_maximo;
    int delta;
    int num_arestas;
    int densidade_arestas;
    int size_clique;
    float delta_clique;
    int num_conflitos;
    int num_cores_busca_local_primeira_melhora_troca_vertice;
    int num_conflitos_busca_local_primeira_melhora_troca_vertice;
    int num_cores_busca_local_melhor_melhora_troca_vertice;
    int num_conflitos_busca_local_melhor_melhora_troca_vertice;
    int num_cores_busca_local_primeira_melhora_cor_menos_frequente;
    int num_conflitos_busca_local_primeira_melhora_cor_menos_frequente;
    int num_cores_busca_local_melhor_melhora_cor_menos_frequente;
    int num_conflitos_busca_local_melhor_melhora_cor_menos_frequente;
    int num_cores_tabu_troca_vertice;
    int num_conflitos_tabu_troca_vertice;
    int num_cores_tabu_menor_frequencia;
    int num_conflitos_tabu_menor_frequencia;
    int num_cores_genetico;
    int num_conflitos_genetico;
};

struct ResponseMaxColor
{
    vector<int> cores;
    int maxColor;
};

struct Solucao
{
    vector<int> solucao;
    int num_conflitos;
    int num_cores;
    double calcularCusto() const
    {
        return (0.9 * num_cores) + (0.1 * num_conflitos);
    }
    bool operator<(const Solucao &other) const
    {
        return calcularCusto() < other.calcularCusto();
    }
    bool operator>(const Solucao &other) const
    {
        return calcularCusto() < other.calcularCusto();
    }
};

struct Selecionados
{
    vector<Solucao> selecionados;
    int num_cores;
    int num_conflitos;
};

struct ResponseBuscaLocal
{
    vector<int> solucao_encontrada;
    int num_conflitos;
    int num_cores;
};

struct ItemGenerateNeighbor
{
    vector<int> vizinho;
    int conflitos;
    int custo;
};

struct ResponseSubstituiCorMenosFrequente
{
    vector<int> solucao;
    vector<int> vertices_modificados;
};

class Problem
{
public:
    vector<int> calcClique(Graph &graph);
    vector<int> troca_cores_aleatorio(vector<int> &cores, Graph &graph);
    vector<int> troca(vector<int> &cores, int pos1, int pos2);
    ResponseSubstituiCorMenosFrequente substitui_cor_menos_frequente(vector<int> &cores, Graph &graph, int cor_menos_frequente, int max_cores);
    ResponseMaxColor getMaxColor(Graph &Graph, int num_vertices, int &grauMaximo);
    ItemGenerateNeighbor tabu_col(Graph &graph, int k, int tabu_size, int rep, int iteracoes, vector<int> &solucao_inicial, int conflito_inicial);
    ItemGenerateNeighbor tabu_col_menor_frequencia(Graph &graph, int k, int tabu_size, int rep, int iteracoes, vector<int> &solucao_inicial, int conflito_inicial);
    HeuristicaConstrutivaResponse heuristicaConstrutiva(Graph &graph, int k);
    int calcConflitoVertice(vector<int> &cores, Graph &graph, int vertice);
    int calcConflitoTotal(vector<int> &cores, Graph &graph);
    int calcCusto(unordered_map<int, int> &frequencias);
    vector<int> calcConflitoGraph(vector<int> &cores, Graph &graph);
    vector<int> ordenaVerticesPorConflito(vector<int> &solucao, Graph &graph);
    unordered_map<int, int> calcFrequenciaCores(vector<int> &solucao);
    vector<int> rankCorPorMenorFrequencia(unordered_map<int, int> &frequencias);
    ResponseBuscaLocal busca_local_melhor_melhora(Graph &graph, vector<int> &solucao_inicial, int iteracoes);
    ResponseBuscaLocal busca_local_primeira_melhora(Graph &graph, vector<int> &solucao_inicial, int iteracoes);
    ResponseBuscaLocal busca_local_primeira_melhora_cor_menos_frequente(Graph &graph, vector<int> &solucao_inicial, int iteracoes);
    ResponseBuscaLocal busca_local_melhor_melhora_cor_menos_frequente(Graph &graph, vector<int> &solucao_inicial, int iteracoes);
    vector<ItemGenerateNeighbor> generate_neighbors_tabu(Graph &graph, vector<int> &current_solution, int k, int rep, const set<pair<int, int>> &tabu_list, int current_conflicts);
    vector<ItemGenerateNeighbor> generate_neighbors_tabu_menor_frequencia(Graph &graph, vector<int> &current_solution, int k, int rep, const set<pair<int, int>> &tabu_list, int current_conflicts);
    void geraPopulacaoInicial(Graph &graph, vector<Solucao> &populacao, int num_solucoes);
    Solucao geraSolucaoAleatoria(Graph &graph);
    Solucao algoritmoGenetico(Graph &graph, vector<Solucao> &populacaoInicial, int num_solucoes, int num_geracoes, int num_mutacoes, int num_selecionados);
    vector<int> torneio(vector<Solucao> &populacao, int num_vencedores);
    Solucao recombinacao(Graph &graph, Solucao &pai1, Solucao &pai2);
    void mutacaoAleatoria(Graph &graph, Solucao &solucao, int max_cor);
};

void Problem::mutacaoAleatoria(Graph &graph, Solucao &solucao, int max_cor)
{
    int vertice = rand() % solucao.solucao.size();
    int cor = (rand() % max_cor);
    solucao.solucao[vertice] = cor;
    solucao.num_conflitos = calcConflitoTotal(solucao.solucao, graph);
    solucao.num_cores = *max_element(solucao.solucao.begin(), solucao.solucao.end());
}

Solucao Problem::recombinacao(Graph &graph, Solucao &pai1, Solucao &pai2)
{
    Solucao filho;
    int ponto_corte_1 = rand() % pai1.solucao.size();
    int ponto_corte_2 = -1;

    while (ponto_corte_2 < ponto_corte_1)
    {
        ponto_corte_2 = rand() % pai1.solucao.size();
    }

    for (int i = 0; i < pai1.solucao.size(); i++)
    {
        if (i < ponto_corte_1 || i > ponto_corte_2)
        {
            filho.solucao.push_back(pai1.solucao[i]);
        }
        else
        {
            filho.solucao.push_back(pai2.solucao[i]);
        }
    }

    filho.num_conflitos = calcConflitoTotal(filho.solucao, graph);
    filho.num_cores = *max_element(filho.solucao.begin(), filho.solucao.end());

    return filho;
}

vector<int> Problem::torneio(vector<Solucao> &populacao, int num_vencedores)
{
    vector<int> vencedores;

    while (vencedores.size() < num_vencedores)
    {
        int index = rand() % populacao.size();
        if (find(vencedores.begin(), vencedores.end(), index) == vencedores.end())
        {
            vencedores.push_back(index);
        }
    }

    return vencedores;
}

Solucao Problem::algoritmoGenetico(Graph &graph, vector<Solucao> &populacaoInicial, int num_solucoes, int num_geracoes, int num_mutacoes, int num_selecionados)
{
    vector<Solucao> populacao = populacaoInicial;
    geraPopulacaoInicial(graph, populacao, num_solucoes);

    for (int i = 0; i < num_geracoes; i++)
    {
        vector<int> selecionados = torneio(populacao, num_selecionados);
        priority_queue<Solucao, vector<Solucao>, greater<Solucao>> heap;

        for (Solucao individuo : populacao)
        {
            heap.push(individuo);
            if (heap.size() > num_solucoes)
            {
                heap.pop();
            }
        }

        for (int selecionado : selecionados)
        {
            for (int j = 0; j < populacao.size(); j++)
            {
                if (selecionado != j)
                {
                    Solucao filho = recombinacao(graph, populacao[selecionado], populacao[j]);
                    heap.push(filho);
                    if (heap.size() > num_solucoes)
                    {
                        heap.pop();
                    }
                }
            }
        }

        vector<Solucao> recombinacoes;
        while (!heap.empty())
        {
            recombinacoes.push_back(heap.top());
            heap.pop();
        }

        vector<bool> mutados(recombinacoes.size(), false);
        for (int j = 0; j < num_mutacoes; j++)
        {
            int index = 0;

            do
            {

                index = rand() % recombinacoes.size();

            } while (mutados[index]);

            mutados[index] = true;
            mutacaoAleatoria(graph, recombinacoes[index], recombinacoes[index].num_cores);
        }

        populacao = recombinacoes;
    }

    return populacao[0];
}

Solucao Problem::geraSolucaoAleatoria(Graph &graph)
{
    srand(time(0));
    int num_vertices = graph.get_num_vertices();
    vector<int> solucao(num_vertices);
    for (int i = 0; i < num_vertices; i++)
    {
        solucao[i] = rand() % (num_vertices);
    }
    Solucao solucao_calculada;
    solucao_calculada.solucao = solucao;
    solucao_calculada.num_conflitos = calcConflitoTotal(solucao, graph);
    solucao_calculada.num_cores = *max_element(solucao.begin(), solucao.end());
    return solucao_calculada;
}

void Problem::geraPopulacaoInicial(Graph &graph, vector<Solucao> &populacao, int num_solucoes)
{
    int qtd_a_gerar = num_solucoes - populacao.size();

    for (int i = 0; i < qtd_a_gerar; i++)
    {
        Solucao solucao = geraSolucaoAleatoria(graph);
        populacao.push_back(solucao);
    }
}
int Problem::calcCusto(unordered_map<int, int> &frequencias)
{
    int custo = 0;
    for (auto &par : frequencias)
    {
        if (par.second > 0)
            custo++;
    }
    return custo;
};

vector<int> Problem::rankCorPorMenorFrequencia(unordered_map<int, int> &frequencias)
{
    unordered_map<int, int> frequencias_aux = frequencias;
    vector<int> rank_cores;

    while (!frequencias_aux.empty())
    {
        auto menorFreqIt = min_element(
            frequencias_aux.begin(), frequencias_aux.end(),
            [](const pair<int, int> &a, const pair<int, int> &b)
            {
                return a.second < b.second;
            });

        if (menorFreqIt->first != 0)
            rank_cores.push_back(menorFreqIt->first);

        frequencias_aux.erase(menorFreqIt);
    }
    return rank_cores;
}

ResponseBuscaLocal Problem::busca_local_melhor_melhora_cor_menos_frequente(Graph &graph, vector<int> &solucao_inicial, int iteracoes)
{
    vector<int> solucao = solucao_inicial;
    vector<int> melhor_vizinho = solucao;

    int melhor_conflito = calcConflitoTotal(solucao, graph);

    vector<int> vertices_ordenados_conflitos = ordenaVerticesPorConflito(solucao, graph);

    unordered_map<int, int> frequenciaCores = calcFrequenciaCores(melhor_vizinho);
    int melhor_custo_global = calcCusto(frequenciaCores);
    int cont_iteracoes = 0;
    while (true)
    {

        vector<int> rank_cores_menor_frequencia = rankCorPorMenorFrequencia(frequenciaCores);
        vector<int> melhor_vizinho_local = melhor_vizinho;
        int melhor_conflito_local = melhor_conflito;
        int melhor_custo_local = melhor_custo_global;
        unordered_map<int, int> frequenciaCores_local = frequenciaCores;
        for (int i = 0; i < rank_cores_menor_frequencia.size(); i++)
        {

            ResponseSubstituiCorMenosFrequente vizinho = substitui_cor_menos_frequente(melhor_vizinho, graph, rank_cores_menor_frequencia[i], melhor_custo_global);
            int conflito_vizinho = calcConflitoTotal(vizinho.solucao, graph);
            unordered_map<int, int> frequenciaVizinho = calcFrequenciaCores(vizinho.solucao);
            int custo_vizinho = calcCusto(frequenciaVizinho);

            if (custo_vizinho < melhor_custo_global && conflito_vizinho <= melhor_conflito)
            {
                melhor_vizinho = vizinho.solucao;
                melhor_conflito = conflito_vizinho;
                melhor_custo_global = custo_vizinho;
                frequenciaCores = frequenciaVizinho;
                melhor_vizinho_local = vizinho.solucao;
                melhor_conflito_local = conflito_vizinho;
                melhor_custo_local = custo_vizinho;
                frequenciaCores_local = frequenciaVizinho;
            }
            else if (conflito_vizinho < melhor_conflito_local)
            {
                melhor_vizinho_local = vizinho.solucao;
                melhor_conflito_local = conflito_vizinho;
                melhor_custo_local = custo_vizinho;
                frequenciaCores_local = frequenciaVizinho;
            }

            cont_iteracoes++;
            if (cont_iteracoes == iteracoes)
            {
                break;
            }
        }

        if (melhor_custo_local <= melhor_custo_global && melhor_conflito_local <= melhor_conflito)
        {
            melhor_vizinho = melhor_vizinho_local;
            melhor_conflito = melhor_conflito_local;
            melhor_custo_global = melhor_custo_local;
            frequenciaCores = frequenciaCores_local;
        }

        if (cont_iteracoes == iteracoes)
        {
            break;
        }
    }

    ResponseBuscaLocal response;
    response.solucao_encontrada = melhor_vizinho;
    response.num_conflitos = melhor_conflito;
    response.num_cores = melhor_custo_global;

    return response;
}

ResponseBuscaLocal Problem::busca_local_primeira_melhora_cor_menos_frequente(Graph &graph, vector<int> &solucao_inicial, int iteracoes)
{
    vector<int> solucao = solucao_inicial;
    vector<int> melhor_vizinho = solucao;

    int melhor_conflito = calcConflitoTotal(solucao, graph);

    vector<int> vertices_ordenados_conflitos = ordenaVerticesPorConflito(solucao, graph);

    unordered_map<int, int> frequenciaCores = calcFrequenciaCores(solucao);

    vector<int> rank_cores_menor_frequencia = rankCorPorMenorFrequencia(frequenciaCores);
    int melhor_custo = calcCusto(frequenciaCores);
    int cont_iteracoes = 0;
    while (cont_iteracoes < iteracoes && melhor_conflito > 0)
    {
        for (int i = 0; i < rank_cores_menor_frequencia.size(); i++)
        {

            ResponseSubstituiCorMenosFrequente vizinho = substitui_cor_menos_frequente(melhor_vizinho, graph, rank_cores_menor_frequencia[i], melhor_custo);
            int conflito_vizinho = calcConflitoTotal(vizinho.solucao, graph);
            unordered_map<int, int> frequenciaVizinho = calcFrequenciaCores(vizinho.solucao);
            int custo_vizinho = calcCusto(frequenciaVizinho);
            if (conflito_vizinho < melhor_conflito || (conflito_vizinho == melhor_conflito && custo_vizinho < melhor_custo))
            {
                melhor_vizinho = vizinho.solucao;
                melhor_conflito = conflito_vizinho;
                melhor_custo = custo_vizinho;
                frequenciaCores = frequenciaVizinho;
                break;
            }
            cont_iteracoes++;
            if (cont_iteracoes == iteracoes)
            {
                break;
            }
        }
    }

    frequenciaCores = calcFrequenciaCores(melhor_vizinho);

    ResponseBuscaLocal response;
    response.solucao_encontrada = melhor_vizinho;
    response.num_conflitos = melhor_conflito;
    response.num_cores = melhor_custo;

    return response;
}

ResponseSubstituiCorMenosFrequente Problem::substitui_cor_menos_frequente(vector<int> &cores, Graph &graph, int cor_menos_frequente, int max_cores)
{
    int num_vertices = graph.get_num_vertices();
    vector<int> cores_aux = cores;
    vector<int> vertices_modificados;

    for (int i = 0; i < num_vertices; ++i)
    {
        if (cores_aux[i] == cor_menos_frequente)
        {
            set<int> cores_vizinhas;
            for (CoupleVertice vizinho : graph.get_neighbors(i))
            {
                cores_vizinhas.insert(cores_aux[vizinho.vertice2]);
            }

            for (int nova_cor = 1; nova_cor <= max_cores; ++nova_cor)
            {
                if (cores_vizinhas.find(nova_cor) == cores_vizinhas.end())
                {
                    cores_aux[i] = nova_cor;
                    vertices_modificados.push_back(i);
                    break;
                }
            }
        }
    }
    ResponseSubstituiCorMenosFrequente response;
    response.solucao = cores_aux;
    response.vertices_modificados = vertices_modificados;
    return response;
};

unordered_map<int, int> Problem::calcFrequenciaCores(vector<int> &solucao)
{
    unordered_map<int, int> frequencias;

    for (int cor : solucao)
    {
        frequencias[cor]++;
    }

    return frequencias;
};

vector<int> Problem::ordenaVerticesPorConflito(vector<int> &solucao, Graph &graph)
{

    vector<int> conflitos = calcConflitoGraph(solucao, graph);

    vector<int> vertices_ordenados_conflitos;

    int maior = INT32_MAX;
    int maior_atual = 0;
    int maior_vertice_atual = 0;
    int pos_ordenacao = 0;

    for (int j = 0; j < graph.get_num_vertices(); j++)
    {
        for (int i = 0; i < conflitos.size(); i++)
        {
            if (conflitos[i] > maior_atual && conflitos[i] < maior)
            {
                maior_atual = conflitos[i];
                maior_vertice_atual = i;
            }
        }

        maior = maior_atual;
        vertices_ordenados_conflitos.push_back(maior_vertice_atual);
    }

    return vertices_ordenados_conflitos;
}

ResponseBuscaLocal Problem::busca_local_primeira_melhora(Graph &graph, vector<int> &solucao_inicial, int iteracoes)
{
    vector<int> solucao = solucao_inicial;
    vector<int> melhor_vizinho = solucao;

    int melhor_conflito = calcConflitoTotal(solucao, graph);
    vector<int> vertices_ordenados_conflitos = ordenaVerticesPorConflito(solucao, graph);

    int cont_iteracoes = 0;
    for (int i = 0; i < vertices_ordenados_conflitos.size(); i++)
    {

        for (CoupleVertice aresta : graph.get_neighbors(i))
        {

            int vertice1 = i;
            int vertice2 = aresta.vertice2;
            vector<int> vizinho = troca(melhor_vizinho, vertice1, vertice2);
            int conflito_vizinho = calcConflitoTotal(vizinho, graph);
            if (conflito_vizinho < melhor_conflito)
            {

                melhor_vizinho = vizinho;
                melhor_conflito = conflito_vizinho;
                break;
            }
            cont_iteracoes++;
            if (cont_iteracoes == iteracoes)
            {
                break;
            }
        }
        if (cont_iteracoes == iteracoes)
        {
            break;
        }
    }

    ResponseBuscaLocal response;
    response.solucao_encontrada = melhor_vizinho;
    response.num_conflitos = melhor_conflito;

    return response;
}

ResponseBuscaLocal Problem::busca_local_melhor_melhora(Graph &graph, vector<int> &solucao_inicial, int iteracoes)
{
    vector<int> solucao = solucao_inicial;
    vector<int> melhor_vizinho = solucao;
    int melhor_custo = *max_element(solucao.begin(), solucao.end());
    int melhor_conflito = calcConflitoTotal(solucao, graph);
    vector<int> vertices_ordenados_conflitos = ordenaVerticesPorConflito(solucao, graph);

    int cont_iteracoes = 0;

    for (int i = 0; i < vertices_ordenados_conflitos.size(); i++)
    {
        vector<int> melhor_vizinho_atual = melhor_vizinho;
        for (CoupleVertice aresta : graph.get_neighbors(i))
        {

            int vertice1 = i;
            int vertice2 = aresta.vertice2;
            vector<int> vizinho = troca(melhor_vizinho_atual, vertice1, vertice2);
            int conflito_vizinho = calcConflitoTotal(vizinho, graph);
            if (conflito_vizinho < melhor_conflito)
            {

                melhor_vizinho = vizinho;
                melhor_conflito = conflito_vizinho;
            }
            cont_iteracoes++;
            if (cont_iteracoes == iteracoes)
            {
                break;
            }
        }

        if (cont_iteracoes == iteracoes)
        {
            break;
        }
    }

    ResponseBuscaLocal response;
    response.solucao_encontrada = melhor_vizinho;
    response.num_conflitos = melhor_conflito;

    return response;
}

vector<ItemGenerateNeighbor> Problem::generate_neighbors_tabu(Graph &graph, vector<int> &solucao_inicial, int k, int rep, const set<pair<int, int>> &tabu_list, int current_conflicts)
{
    vector<int> solucao = solucao_inicial;
    vector<int> melhor_vizinho = solucao;
    int melhor_custo = k;
    int melhor_conflito = current_conflicts;
    vector<int> vertices_ordenados = ordenaVerticesPorConflito(solucao, graph);
    vector<ItemGenerateNeighbor> vizinhos;

    for (int vertice : vertices_ordenados)
    {
        if (vizinhos.size() >= rep)
            break;

        for (CoupleVertice aresta : graph.get_neighbors(vertice))
        {
            if (vizinhos.size() >= rep)
                break;

            vector<int> vizinho = troca(solucao, vertice, aresta.vertice2);

            int new_conflicts = calcConflitoTotal(vizinho, graph);

            if (tabu_list.find({vertice, vizinho[vertice]}) == tabu_list.end() || new_conflicts < current_conflicts)
            {
                ItemGenerateNeighbor vizinho_item;
                vizinho_item.vizinho = vizinho;
                vizinho_item.conflitos = new_conflicts;
                vizinhos.push_back(vizinho_item);
            }
        }
    }

    return vizinhos;
}

vector<ItemGenerateNeighbor> Problem::generate_neighbors_tabu_menor_frequencia(Graph &graph, vector<int> &solucao_inicial, int k, int rep, const set<pair<int, int>> &tabu_list, int current_conflicts)
{
    vector<int> solucao = solucao_inicial;
    vector<int> melhor_vizinho = solucao;
    int melhor_custo = k;
    int melhor_conflito = current_conflicts;
    vector<int> vertices_ordenados = ordenaVerticesPorConflito(solucao, graph);
    unordered_map<int, int> frequenciaCores = calcFrequenciaCores(melhor_vizinho);

    vector<int> rank_cores_menor_frequencia = rankCorPorMenorFrequencia(frequenciaCores);
    vector<ItemGenerateNeighbor> vizinhos;

    for (int i = 0; i < rank_cores_menor_frequencia.size(); i++)
    {
        if (vizinhos.size() >= rep)
            break;

        ResponseSubstituiCorMenosFrequente vizinho = substitui_cor_menos_frequente(solucao, graph, rank_cores_menor_frequencia[i], melhor_custo);
        int conflito_vizinho = calcConflitoTotal(vizinho.solucao, graph);
        unordered_map<int, int> frequenciaVizinho = calcFrequenciaCores(vizinho.solucao);
        int custo_vizinho = calcCusto(frequenciaVizinho);

        bool tabu = false;
        for (int vertice : vizinho.vertices_modificados)
        {
            if (tabu_list.find({vertice, vizinho.solucao[vertice]}) != tabu_list.end())
            {
                tabu = true;
                break;
            }
        }

        if (!tabu || conflito_vizinho < current_conflicts || (custo_vizinho < melhor_custo && conflito_vizinho < current_conflicts))
        {
            ItemGenerateNeighbor vizinho_item;
            vizinho_item.vizinho = vizinho.solucao;
            vizinho_item.conflitos = conflito_vizinho;
            vizinho_item.custo = custo_vizinho;
            vizinhos.push_back(vizinho_item);
        }
    }

    return vizinhos;
}

ItemGenerateNeighbor Problem::tabu_col(Graph &graph, int k, int tabu_size, int rep, int iteracoes, vector<int> &solucao_inicial, int conflito_inicial)
{
    srand(time(0));

    ItemGenerateNeighbor vizinho_atual;
    vizinho_atual.vizinho = solucao_inicial;
    vizinho_atual.conflitos = conflito_inicial;
    vizinho_atual.custo = k;

    deque<pair<int, int>> tabu_list;
    int nbiter = 0;

    while (vizinho_atual.conflitos > 0 && nbiter < iteracoes)
    {
        vector<ItemGenerateNeighbor> neighbors = generate_neighbors_tabu(graph, vizinho_atual.vizinho, k, rep, {tabu_list.begin(), tabu_list.end()}, vizinho_atual.conflitos);

        ItemGenerateNeighbor melhor_vizinho;
        melhor_vizinho.vizinho = vizinho_atual.vizinho;
        melhor_vizinho.conflitos = numeric_limits<int>::max();

        for (const ItemGenerateNeighbor &neighbor : neighbors)
        {
            if (neighbor.conflitos < melhor_vizinho.conflitos)
            {
                melhor_vizinho = neighbor;
            }
        }

        if (vizinho_atual.vizinho != melhor_vizinho.vizinho)
        {
            for (int i = 0; i < graph.get_num_vertices(); ++i)
            {
                if (vizinho_atual.vizinho[i] != melhor_vizinho.vizinho[i])
                {
                    tabu_list.push_back({i, melhor_vizinho.vizinho[i]});
                    if (tabu_list.size() > tabu_size)
                    {
                        tabu_list.pop_front();
                    }
                    break;
                }
            }
        }

        vizinho_atual = melhor_vizinho;

        nbiter++;
    }

    return vizinho_atual;
}

ItemGenerateNeighbor Problem::tabu_col_menor_frequencia(Graph &graph, int k, int tabu_size, int rep, int iteracoes, vector<int> &solucao_inicial, int conflito_inicial)
{
    srand(time(0));

    ItemGenerateNeighbor vizinho_atual;
    vizinho_atual.vizinho = solucao_inicial;
    vizinho_atual.conflitos = conflito_inicial;
    vizinho_atual.custo = k;

    deque<pair<int, int>> tabu_list;
    int nbiter = 0;

    while (nbiter < iteracoes)
    {
        vector<ItemGenerateNeighbor> neighbors = generate_neighbors_tabu_menor_frequencia(graph, vizinho_atual.vizinho, vizinho_atual.custo, rep, {tabu_list.begin(), tabu_list.end()}, vizinho_atual.conflitos);

        ItemGenerateNeighbor melhor_vizinho;
        melhor_vizinho.vizinho = vizinho_atual.vizinho;
        melhor_vizinho.conflitos = numeric_limits<int>::max();

        for (const ItemGenerateNeighbor &neighbor : neighbors)
        {
            if (neighbor.custo < melhor_vizinho.custo && neighbor.conflitos < melhor_vizinho.conflitos)
            {
                melhor_vizinho = neighbor;
            }
            else if (neighbor.conflitos < melhor_vizinho.conflitos)
            {
                melhor_vizinho = neighbor;
            }
        }

        if (vizinho_atual.vizinho != melhor_vizinho.vizinho)
        {
            for (int i = 0; i < graph.get_num_vertices(); ++i)
            {
                if (vizinho_atual.vizinho[i] != melhor_vizinho.vizinho[i])
                {
                    tabu_list.push_back({i, melhor_vizinho.vizinho[i]});
                    if (tabu_list.size() > tabu_size)
                    {
                        tabu_list.pop_front();
                    }
                    break;
                }
            }
        }

        vizinho_atual = melhor_vizinho;

        nbiter++;
    }

    return vizinho_atual;
}

int Problem::calcConflitoTotal(vector<int> &cores, Graph &graph)
{
    int conflitos = 0;

    for (int vertice = 0; vertice < graph.get_num_vertices(); vertice++)
    {
        conflitos += calcConflitoVertice(cores, graph, vertice);
    }
    return conflitos / 2;
}

vector<int> Problem::calcConflitoGraph(vector<int> &cores, Graph &graph)
{
    vector<int> conflitos(graph.get_num_vertices(), 0);

    for (int vertice = 0; vertice < graph.get_num_vertices(); vertice++)
    {
        conflitos[vertice] = calcConflitoVertice(cores, graph, vertice);
    }
    return conflitos;
}

int Problem::calcConflitoVertice(vector<int> &cores, Graph &graph, int vertice)
{
    int conflito = 0;
    for (CoupleVertice vizinho : graph.get_neighbors(vertice))
    {
        if (cores[vizinho.vertice2] == cores[vertice])
            conflito++;
    }
    return conflito;
}

vector<int> Problem::troca_cores_aleatorio(vector<int> &cores, Graph &graph)
{
    vector<int> cores_aux = cores;
    srand(time(0));
    int vertice1 = rand() % graph.get_num_vertices();
    int vertice2 = 0;
    do
    {
        vertice2 = rand() % graph.get_num_vertices();
    } while (vertice1 == vertice2);

    swap(cores_aux[vertice1], cores_aux[vertice2]);

    return cores_aux;
}

vector<int> Problem::troca(vector<int> &cores, int pos1, int pos2)
{
    vector<int> cores_aux = cores;
    swap(cores_aux[pos1], cores_aux[pos2]);
    return cores_aux;
}

vector<int> Problem::calcClique(Graph &graph)
{
    int n = graph.get_num_vertices();
    vector<int> vertices(n);
    for (int i = 0; i < n; ++i)
        vertices[i] = i;

    sort(vertices.begin(), vertices.end(), [&](int a, int b)
         { return graph.get_neighbors(a).size() > graph.get_neighbors(b).size(); });

    vector<int> clique;
    for (int v : vertices)
    {
        bool canAdd = true;
        for (int u : clique)
        {
            if (!graph.ehVizinho(u, v))
            {
                canAdd = false;
                break;
            }
        }
        if (canAdd)
            clique.push_back(v);
    }
    return clique;
}

ResponseMaxColor Problem::getMaxColor(Graph &graph, int num_vertices, int &grauMaximo)
{
    vector<int> cores(num_vertices, -1);
    vector<int> saturacao(num_vertices, 0);
    vector<int> graus(num_vertices, 0);

    for (int i = 0; i < num_vertices; i++)
    {
        graus[i] = graph.get_neighbors(i).size();
        if (graus[i] > grauMaximo)
        {
            grauMaximo = graus[i];
        }
    }

    for (int passada = 0; passada < num_vertices; passada++)
    {
        int melhorVertice = -1;
        int maiorSaturacao = -1;
        int maiorGrau = -1;

        for (int vertice = 0; vertice < num_vertices; vertice++)
        {
            if (cores[vertice] == -1)
            {
                if (saturacao[vertice] > maiorSaturacao ||
                    (saturacao[vertice] == maiorSaturacao && graus[vertice] > maiorGrau))
                {
                    melhorVertice = vertice;
                    maiorSaturacao = saturacao[vertice];
                    maiorGrau = graus[vertice];
                }
            }
        }

        int menorCor = 1;
        set<int> coresIndisponiveis;
        for (CoupleVertice vizinho : graph.get_neighbors(melhorVertice))
        {
            if (cores[vizinho.vertice2] != -1)
            {
                coresIndisponiveis.insert(cores[vizinho.vertice2]);
            }
        }

        while (coresIndisponiveis.count(menorCor) > 0)
        {
            menorCor++;
        }

        cores[melhorVertice] = menorCor;
        for (CoupleVertice vizinho : graph.get_neighbors(melhorVertice))
        {
            saturacao[vizinho.vertice2]++;
        }
    }

    int maxCor = 0;
    for (int i = 0; i < num_vertices; i++)
    {
        if (cores[i] > maxCor)
        {
            maxCor = cores[i];
        }
    }

    ResponseMaxColor response;
    response.cores = cores;
    response.maxColor = maxCor;

    return response;
}

HeuristicaConstrutivaResponse Problem::heuristicaConstrutiva(Graph &graph, int initialK = 10)
{
    HeuristicaConstrutivaResponse response;
    vector<int> clique = calcClique(graph);
    int k = clique.size();

    int num_arestas = 0;
    for (int i = 0; i < graph.get_num_vertices(); i++)
    {
        num_arestas += graph.get_neighbors(i).size();
    }

    int num_vertices = graph.get_num_vertices();
    int grauMaximo = 0;
    cout << "executando construtiva" << endl;
    ResponseMaxColor responseMaxColor = getMaxColor(graph, num_vertices, grauMaximo);
    int maxCor = responseMaxColor.maxColor;

    int conflitos = calcConflitoTotal(responseMaxColor.cores, graph);

    response.num_conflitos = conflitos;
    vector<Solucao> populacaoInicialSolucoes;

    Solucao solucao_heuristica;
    solucao_heuristica.solucao = responseMaxColor.cores;
    solucao_heuristica.num_conflitos = conflitos;
    solucao_heuristica.num_cores = maxCor;
    populacaoInicialSolucoes.push_back(solucao_heuristica);
    //////////////////////////////////////////////////////////////////////////////////////////////
    // ItemGenerateNeighbor resposta_tabu = tabu_col(graph, maxCor, 30, 10, 10000, responseMaxColor.cores, conflitos);
    // response.num_conflitos_tabu_troca_vertice = resposta_tabu.conflitos;
    // response.num_cores_tabu_troca_vertice = resposta_tabu.custo;

    // Solucao solucao_tabu;
    // solucao_tabu.solucao = resposta_tabu.vizinho;
    // solucao_tabu.num_conflitos = resposta_tabu.conflitos;
    // solucao_tabu.num_cores = resposta_tabu.custo;
    // populacaoInicalSolucoes.push_back(solucao_tabu);
    //////////////////////////////////////////////////////////////////////////////////////////////
    // resposta_tabu = tabu_col_menor_frequencia(graph, maxCor, 30, 10, 10000, responseMaxColor.cores, conflitos);
    // response.num_conflitos_tabu_menor_frequencia = resposta_tabu.conflitos;
    // response.num_cores_tabu_menor_frequencia = resposta_tabu.custo;

    // Solucao solucao_tabu_menor_frequencia;
    // solucao_tabu_menor_frequencia.solucao = resposta_tabu.vizinho;
    // solucao_tabu_menor_frequencia.num_conflitos = resposta_tabu.conflitos;
    // solucao_tabu_menor_frequencia.num_cores = resposta_tabu.custo;
    // // populacaoInicialSolucoes.push_back(solucacao_tabu_menor_frequencia);
    // //////////////////////////////////////////////////////////////////////////////////////////////
    // ResponseBuscaLocal busca = busca_local_primeira_melhora(graph, responseMaxColor.cores, 10000);

    // response.num_conflitos_busca_local_primeira_melhora_troca_vertice = busca.num_conflitos;
    // response.num_cores_busca_local_primeira_melhora_troca_vertice = maxCor;

    // Solucao solucao_busca_local_primeira_melhora;
    // solucao_busca_local_primeira_melhora.solucao = busca.solucao_encontrada;
    // solucao_busca_local_primeira_melhora.num_conflitos = busca.num_conflitos;
    // solucao_busca_local_primeira_melhora.num_cores = busca.num_cores;
    // // populacaoInicialSolucoes.push_back(solucao_busca_local_primeira_melhora);

    // //////////////////////////////////////////////////////////////////////////////////////////////
    // busca = busca_local_primeira_melhora_cor_menos_frequente(graph, responseMaxColor.cores, 20000);

    // response.num_conflitos_busca_local_primeira_melhora_cor_menos_frequente = busca.num_conflitos;
    // response.num_cores_busca_local_primeira_melhora_cor_menos_frequente = busca.num_cores;

    // Solucao solucao_busca_local_primeira_melhora_cor_menos_frequente;
    // solucao_busca_local_primeira_melhora_cor_menos_frequente.solucao = busca.solucao_encontrada;
    // solucao_busca_local_primeira_melhora_cor_menos_frequente.num_conflitos = busca.num_conflitos;
    // solucao_busca_local_primeira_melhora_cor_menos_frequente.num_cores = busca.num_cores;
    // populacaoInicialSolucoes.push_back(solucao_busca_local_primeira_melhora_cor_menos_frequente);

    // int n_cor_busca_pri_m_c_f = busca.num_cores;

    //////////////////////////////////////////////////////////////////////////////////////////////

    // busca = busca_local_melhor_melhora_cor_menos_frequente(graph, responseMaxColor.cores, 20000);
    // response.num_conflitos_busca_local_melhor_melhora_cor_menos_frequente = busca.num_conflitos;
    // response.num_cores_busca_local_melhor_melhora_cor_menos_frequente = busca.num_cores;

    // Solucao solucao_busca_local_melhor_melhora_cor_menos_frequente;
    // solucao_busca_local_melhora_cor_menos_frequente.solucao = busca.solucao_encontrada;
    // solucao_busca_local_melhora_cor_menos_frequente.num_conflitos = busca.num_conflitos;
    // solucao_busca_local_melhora_cor_menos_frequente.num_cores = busca.num_cores;
    // populacaoInicialSolucoes.push_back(solucao_busca_local_melhor_melhora_cor_menos_frequente);

    // int n_cor_busca_melhor_m_c_m_f = busca.num_cores;

    //////////////////////////////////////////////////////////////////////////////////////////////

    // busca = busca_local_melhor_melhora(graph, responseMaxColor.cores, 10000);

    // response.num_conflitos_busca_local_melhor_melhora_troca_vertice = busca.num_conflitos;
    // response.num_cores_busca_local_melhor_melhora_troca_vertice = maxCor;
    // Solucao solucao_busca_local_melhor_melhora;
    // solucao_busca_local_melhor_melhora.solucao = busca.solucao_encontrada;
    // solucao_busca_local_melhor_melhora.num_conflitos = busca.num_conflitos;
    // solucao_busca_local_melhor_melhora.num_cores = busca.num_cores;
    // populacaoInicialSolucoes.push_back(solucao_busca_local_melhor_melhora);

    //////////////////////////////////////////////////////////////////////////////////////////////

    int num_solucoes = 20;
    int num_geracoes = 10000;
    int num_mutacoes = 4;
    int num_selecionados = 6;

    cout << "executando genético" << endl;
    auto inicio = high_resolution_clock::now();
    Solucao solucaoGenetica = algoritmoGenetico(graph, populacaoInicialSolucoes, num_solucoes, num_geracoes, num_mutacoes, num_selecionados);
    response.num_cores_genetico = solucaoGenetica.num_cores;
    response.num_conflitos_genetico = solucaoGenetica.num_conflitos;
    auto fim = high_resolution_clock::now();

    // Calcula a duração
    auto duracao = duration_cast<milliseconds>(fim - inicio);

    cout << "Tempo de execução: " << duracao.count() << " ms" << endl;
    // cout << "Solucao Genetico: " << solucaoGenetica.num_cores << " " << solucaoGenetica.num_conflitos << endl;

    response.size_clique = k;
    response.delta_clique = (float)(maxCor - k) / k;
    response.num_cores_encontradas = maxCor;
    response.num_grau_maximo = grauMaximo;
    response.delta = grauMaximo - maxCor;
    response.num_arestas = num_arestas;
    response.densidade_arestas = (float)(num_arestas) / num_vertices;
    response.num_vertices = num_vertices;

    return response;
}
#endif
