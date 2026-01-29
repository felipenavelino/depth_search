import pickle
import numpy as np
import plotly.graph_objects as go
import matplotlib.pyplot as plt
import os


def busca_profundidade_rotular(matriz, x, y, z, visitado, direcoes, rotulo, valor_celula):
    pilha = [(x, y, z)]
    visitado[x][y][z] = True
    matriz[x][y][z] = rotulo
    tamanho = 1
    voxels_componente = [(x, y, z)]

    while pilha:
        cx, cy, cz = pilha.pop()
        for dx, dy, dz in direcoes:
            nx, ny, nz = cx + dx, cy + dy, cz + dz
            if 0 <= nx < len(matriz) and 0 <= ny < len(matriz[0]) and 0 <= nz < len(matriz[0][0]):
                if not visitado[nx][ny][nz] and matriz[nx][ny][nz] == valor_celula:
                    visitado[nx][ny][nz] = True
                    matriz[nx][ny][nz] = rotulo
                    pilha.append((nx, ny, nz))
                    tamanho += 1
                    voxels_componente.append((nx, ny, nz))

    return tamanho, voxels_componente


def escolher_vizinhanca():
    print("Escolha o tipo de vizinhança para a análise:")
    print("1. Vizinhança-6")
    print("2. Vizinhança-18")
    print("3. Vizinhança-26")
    
    while True:
        escolha = input("Digite o número correspondente à vizinhança (1/2/3): ")
        if escolha == '1':
            return [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)]
        elif escolha == '2':
            direcoes = [(dx, dy, dz) for dx in [-1, 0, 1] for dy in [-1, 0, 1] for dz in [-1, 0, 1] if (dx, dy, dz) != (0, 0, 0)]
            return [d for d in direcoes if sum(abs(i) for i in d) <= 2]
        elif escolha == '3':
            return [(dx, dy, dz) for dx in [-1, 0, 1] for dy in [-1, 0, 1] for dz in [-1, 0, 1] if (dx, dy, dz) != (0, 0, 0)]
        else:
            print("Opção inválida. Tente novamente.")


def rotular_componentes_conectados(matriz, direcoes):
    visitado = [[[False for _ in range(len(matriz[0][0]))] for _ in range(len(matriz[0]))] for _ in range(len(matriz))]
    rotulo = 2
    tipos_celulas = [(255, "proliferativa"), (200, "quiescente"), (140, "necrótica")]
    tamanhos_componentes = {"proliferativa": [], "quiescente": [], "necrótica": []}

    for valor_celula, nome_celula in tipos_celulas:
        print(f"Rotulando componentes conectados para células {nome_celula} (valor {valor_celula})...")
        for x in range(len(matriz)):
            for y in range(len(matriz[0])):
                for z in range(len(matriz[0][0])):
                    if not visitado[x][y][z] and matriz[x][y][z] == valor_celula:
                        tamanho, voxels = busca_profundidade_rotular(matriz, x, y, z, visitado, direcoes, rotulo, valor_celula)
                        rotulo += 1
                        tamanhos_componentes[nome_celula].append((tamanho, voxels))

    return matriz, tamanhos_componentes


def contar_tipos_celulas(tamanhos_componentes):
    contagem_celulas = {}
    for tipo_celula, componentes in tamanhos_componentes.items():
        total = sum(tamanho for tamanho, _ in componentes)
        contagem_celulas[tipo_celula] = total
    return contagem_celulas


def salvar_histograma(contagem_celulas, diretorio_saida="saida_histograma"):
    os.makedirs(diretorio_saida, exist_ok=True)
    
    plt.figure(figsize=(10, 6))
    barras = plt.bar(contagem_celulas.keys(), contagem_celulas.values(), color='skyblue', edgecolor='black')
    plt.xlabel('Tipos de Células', fontsize=14)
    plt.ylabel('Quantidade', fontsize=14)
    plt.title('Histograma de Quantidade de Células por Tipo', fontsize=16)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.grid(axis='y', linestyle='--', alpha=0.7)

    for barra in barras:
        valor_y = barra.get_height()
        plt.text(barra.get_x() + barra.get_width()/2, valor_y + 0.5, int(valor_y), ha='center', va='bottom', fontsize=12)

    caminho_saida = os.path.join(diretorio_saida, "histograma_celulas_preciso.png")
    plt.savefig(caminho_saida, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Histograma com precisão aprimorada salvo em: {caminho_saida}")


def plotar_maior_componente(tamanhos_componentes, diretorio_saida="visualizacoes_3D"):
    os.makedirs(diretorio_saida, exist_ok=True)
    for tipo_celula, componentes in tamanhos_componentes.items():
        if componentes:
            maior_componente = max(componentes, key=lambda x: x[0])
            tamanho, voxels = maior_componente
            x, y, z = zip(*voxels)
            fig = go.Figure(data=[go.Scatter3d(
                x=x, y=y, z=z,
                mode='markers',
                marker=dict(
                    size=4,
                    color=np.full(len(x), 1),
                    opacity=0.8
                )
            )])
            fig.update_layout(
                title=f"Maior Componente de Células {tipo_celula}",
                scene=dict(
                    xaxis_title='X',
                    yaxis_title='Y',
                    zaxis_title='Z'
                )
            )
            caminho_saida = os.path.join(diretorio_saida, f"maior_componente_{tipo_celula}.png")
            fig.write_image(caminho_saida)
            print(f"Visualização do maior componente de células {tipo_celula} salva em: {caminho_saida}")


def plotar_componentes_3D(matriz, tipos_celulas, diretorio_saida="visualizacoes_3D"):
    os.makedirs(diretorio_saida, exist_ok=True)
    x, y, z = np.nonzero(matriz)
    valores = matriz[x, y, z]
    valores_validos = [v for v, _ in tipos_celulas]
    mascara = np.isin(valores, valores_validos)
    x, y, z, valores = x[mascara], y[mascara], z[mascara], valores[mascara]
    fig = go.Figure(data=[go.Scatter3d(
        x=x, y=y, z=z,
        mode='markers',
        marker=dict(
            size=3,
            color=valores,
            colorscale='Jet',
            opacity=0.8,
            colorbar=dict(title="Rótulo")
        )
    )])
    fig.update_layout(scene=dict(
        xaxis_title='X',
        yaxis_title='Y',
        zaxis_title='Z'
    ))
    caminho_saida = os.path.join(diretorio_saida, "componentes_3D.png")
    fig.write_image(caminho_saida)
    print(f"Visualização da matriz 3D rotulada salva em: {caminho_saida}")


def salvar_distribuicao_componentes(tamanhos_componentes, arquivo_saida="distribuicao_componentes.txt"):
    with open(arquivo_saida, 'w') as f:
        for tipo_celula, componentes in tamanhos_componentes.items():
            tamanhos = [tamanho for tamanho, _ in componentes]
            f.write(f"Células {tipo_celula}:\n")
            f.write(f"Número de agrupamentos: {len(tamanhos)}\n")
            f.write("Tamanhos dos agrupamentos: " + ", ".join(map(str, sorted(tamanhos))) + "\n\n")
    print(f"Distribuição dos tamanhos dos componentes salva em: {arquivo_saida}")


with open('matrix_3d', 'rb') as f:
    matriz_3d = pickle.load(f)


direcoes = escolher_vizinhanca()


matriz_rotulada, tamanhos_componentes = rotular_componentes_conectados(matriz_3d, direcoes)


contagem_celulas = contar_tipos_celulas(tamanhos_componentes)


salvar_histograma(contagem_celulas)


salvar_distribuicao_componentes(tamanhos_componentes)


plotar_maior_componente(tamanhos_componentes)


tipos_celulas = [(255, "proliferativa"), (200, "quiescente"), (140, "necrótica")]
plotar_componentes_3D(matriz_rotulada, tipos_celulas)

