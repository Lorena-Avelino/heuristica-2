import random

def gerar_grafo_aleatorio(num_vertices):

    grafo = {}
    for i in range(num_vertices):
        grafo[i] = []

    num_arestas_possiveis = num_vertices * (num_vertices - 1) // 2
    num_arestas_alvo = int(0.5 * num_arestas_possiveis)

    arestas_geradas = 0
    for i in range(num_vertices):
        for j in range(i + 1, num_vertices):
            if random.random() < 0.5 and arestas_geradas < num_arestas_alvo:  # Probabilidade de 0.5 para cada aresta
                grafo[i].append(j)
                grafo[j].append(i)
                arestas_geradas += 1

    return grafo

def salvar_grafo_txt(config, nome_arquivo):
    with open(nome_arquivo, 'w') as arquivo:
        id_graph = 0
        for type_graph in config:
            num_graphs = type_graph[0]
            num_vertices = type_graph[1]

            for i in range(num_graphs):
                grafo = gerar_grafo_aleatorio(num_vertices)
                for vertice, vizinhos in grafo.items():
                    linha = f"{id_graph}: {vertice}: {', '.join(map(str, vizinhos))},\n"
                    arquivo.write(linha)
                id_graph += 1

# Exemplo de uso

config_random_graphs = [
[20,100],
[10,300],
[5,500],
[2,1000]
];

nome_arquivo = "grafo_aleatorio.txt"
salvar_grafo_txt(config_random_graphs, nome_arquivo)