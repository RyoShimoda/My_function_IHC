# My_function_IHC
## はじめに

このリポジトリには、Immunohistochemistryのデータ分析とグラフの描画を効率化するためのR自作関数が格納されています。コードを実行すると、自動的にグラフが描画・保存され、解析用のデータセット (csv) が出力されます。Immunohistochemistryにおける切片の写真には、 (BZ-X710, Keyence) を用い、解析にはBZ-X analyzer software (Keyence) を用いることを想定しています。

## 概要

この関数は、Immunohistochemistryの膨大なデータ量を解析するのに大幅な時間がかかる課題を解決するために作成しました。

- 特徴1 データの入ったフォルダを選択するだけで、すべてのデータを読み込み、整形したデータを出力します。
- 特徴2 ggplot2ベースのグラフが細胞カウント (Count)、面積 (Area)、輝度 (Luminance)、(赤) 色輝度 (Red Intensity) でそれぞれ描画されるだけでなく、検体ごとの個別データ、ヒストグラムも出力されます。
- 特徴3 引数で自身の好みの描画幅や高さ、線の太さや色等も変更できます。

## 使い方
### データセットの用意
この関数は前頭前皮質 (Pre-Frontal Cortex: PFC) と海馬 (Hippocampus) を二色（活性細胞を赤、NeuNを緑）で解析することを想定しています。各組織の40μm切片をBZ-X710で6枚撮影し、データファイルに群の名前 (安静:SED等)、領域 (PL等)、撮影位置 (1-6等)、色 (R、G) を入れます (例: MacroSED1_PL_1_R.csv)。海馬は9領域 (dDGsp, dDGip, dCA1, dCA2, dCA3, vDGsp, vDGip, vCA1, vCA3) で、背側 (dorsal) のCA2以外は各4枚、dCA2と腹側 (ventral) は各2枚ずつ解析することを想定しています。また、実施した実験の海馬の活性細胞数が少ない場合、この関数で出力される図はあまり参考にならないため、活性細胞数を手動カウントしたデータを用意し、この関数で出力される整形データと合わせて解析を行う場合もあります (Hiraga, Shimoda _et al_ 2026)。BZ-X analyzer softwareで緑光40、赤光20等で出力した実験データを、以下のようにフォルダ (G40_R20等) に入れ、それをPFC等のフォルダに格納してください。

```
└── PFC/            # このフォルダにResultsフォルダが生成され、結果が格納される。
    └── G40_R20/    # 解析するデータファイルが入ったフォルダ。
 ```

### 関数の読み込み

GitHubからスクリプトファイルを読み込むことで関数が利用可能になります。 必要なライブラリーも関数の中で読み込まれます。

```R
# GitHubから直接スクリプトを読み込む
source("https://github.com/RyoShimoda/My_function_IHC/blob/main/MF_IHC.R")
```

## 使用例 

https://github.com/user-attachments/assets/c3727a3b-9d31-4272-a051-3b4cc4a378ff

この関数を使用すると、Folder_nameで指定したフォルダにResultsのフォルダが生成され、結果が格納されます。
```
└── PFC/                   # データセットが入ったフォルダ、引数folder_nameで指定。作業ディレクトリがここの場合は指定する必要はない。
    ├── G40_R20/           # 解析するデータが格納されているフォルダ
    │   ├── ○○.csv
    │   └── ⋮
    └── Results/           # 解析結果が格納されるフォルダ
        ├── Plots/         # 描画されたグラフが格納されるフォルダ
        │   ├── All/       # すべての領域の棒グラフが横並びになるグラフが保存
        │   │   ├── ○○.png # 細胞数、面積、輝度、赤色輝度がそれぞれ保存
        │   │   └── ⋮
        │   ├── DV/        # 背側、腹側に分かれたグラフが出力。PFCの場合は出力されない
        │   ├── Histogram/ # 最大輝度と平均輝度のヒストグラムが保存
        │   │   ├── ○○.png
        │   │   └── ⋮
        │   ├── Individual_data/ # 細胞数、面積、輝度、赤色輝度それぞれの個別データが保存
        │   │   ├── ○○.png
        │   │   └── ⋮
        │   ├── Region/    # 領域ごとの細胞数、面積、輝度、赤色輝度の棒グラフが保存
        │   │   ├── ○○.png
        │   │   └── ⋮
        ├── GatheringData.csv    # 整形したデータが入ったcsvファイル
        └── HistData.csv         # ヒストグラム用に整形したデータが入ったcsvファイル
```

## 関数の詳細

| 引数名 | デフォルト値 | 説明 |
| :--- | :---: | :--- |
| `Folder_name` | "" | ここで指定したフォルダにResultsフォルダを生成し、グラフを保存 |
| `Region_list`| (必須) | 領域名のリスト。自身でリストを作り指定 (例: PFCList <- c("PL", "IL")) |
| `Path` | (必須) | データファイルが入っているフォルダまでのpathを指定 |
| `Red_name` | "RED" | 赤色の細胞の名前 |
| `Green_name` | "GREEN" | 緑色の細胞の名前 |
| `PFC` | "F" | PFCの場合は PFC = "" のように指定 |
| `D_V` | "T" | PFCの場合は D_V = "" のように指定 |
| `number_of_region` | 9 | PFCの場合は number_of_region = 2 のように指定 |
&#8942;

よく使う引数は上記の通りです。この関数はデフォルトで海馬解析用の引数が指定されています。使用例の動画はPFCのデータを解析しています。データセットのフォルダに作業ディレクトリを指定していないため、Folder_name = "Demo"と指定し、赤色がc-Fos、緑色がNeuNを表しているので、Red_name = "c-Fos"、Green_name = "NeuN"と指定し、グラフに描画される領域の数が二つと少ないので、それぞれ引数でグラフの縦横の幅を指定しています。<br>
海馬の解析では以下のような最短コードで実行できます。引数ではその他、群の名前、グラフの各線やテキストの大きさ、群の色などは個別に指定できます。
```R
# データセットが入っているフォルダを作業ディレクトリにしている場合はfolder_nameも省略可能
# folder_nameで指定したフォルダの下に、グラフ等が格納されるResultsフォルダが生成されます
# 赤色、緑色細胞の名前がデフォルトのまま出力されます
Region_Hippocampus <- c("dDGsp", "dDGip", "dCA3", "dCA2", "dCA1", "vDGsp", "vDGip", "vCA3", "vCA1")
mf_IHC(Region_list = Region_Hippocampus, Path = "パス")
```
群の名前は (Shimoda _et al_., 2024) のとおり安静群はSED、低強度運動群はLIEとしているので、引数のデフォルトは統制群がSED = "SED"、処理群がLIE  ="LIE"として指定されていますが、群の名前を変える場合は統制群をSED  ="Control"、処理群をLIE = "EXE"等と変更してください。その場合、graphic_reference = "(統制群の名前)"の引数も追加してください。ただ、上記の使用例の箇所でも言及しているとおり、データファイルの名前にその群の名前を入れてください。また、群ごとの色分けも引数で変更可能です。

## 活用実績
この関数は以下の論文、リポジトリでのデータ解析に活用されています。
- Article: [Regular light-intensity exercise accelerates contextual fear extinction with reduced dorsal CA3 activation in male rats. Neurochemistry International,Volume 198,2026,106191](https://www.sciencedirect.com/science/article/pii/S0197018626000823?via%3Dihub)<br>
Repository: [Analyse_Neurochem.-Int.](https://github.com/RyoShimoda/Analyse_Neurochem.-Int.)

## 動作環境

* R バージョン 4.6.1
* 使用パッケージ: ggplot2, dplyr, tidyr, stringr, patchwork
