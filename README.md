# kits

kits 把不同年齡、答不同題目的孩子放到同一把尺上，並做出同齡常模。
它是 KIT（Kids in Taiwan，台灣幼兒發展調查）發展量尺的估計流程：校準、世代連結、接到 GSED 的 D 分數、計分、常模、檢查。
模型只寫一次，導數都由 RTMB 自動微分（校準 CML、計分 WLE、門檻的 SE、常模 shash）；封閉解用 base R；整理資料用 tidyverse 的 dplyr、tidyr、purrr。

第一次用：照下面裝好，跑「第一個例子」，再讀入門篇 `vignette("kits-flow", package = "kits")`。

## 安裝

先裝需要的套件：

```r
install.packages(c("RTMB", "dplyr", "tidyr", "purrr", "plotly", "htmlwidgets"))
```

RTMB 在 Linux 上要編譯，第一次可能要好幾分鐘。

再裝 kits，從 GitHub：

```r
remotes::install_github("jiewenTsai/kits")
```

或用下載的原始碼套件：

```r
install.packages("kits_0.25.1.tar.gz", repos = NULL, type = "source")
```

不能寫系統的 library 時（例如學校的電腦），裝到自己的資料夾：

```r
dir.create("~/R/mylib", recursive = TRUE)
install.packages("kits_0.25.1.tar.gz", repos = NULL, type = "source", lib = "~/R/mylib")
.libPaths(c("~/R/mylib", .libPaths()))
library(kits)
```

tar.gz 裡已經有建好的教學。
從原始碼的資料夾安裝、又沒有建好的教學時，才需要先 `R CMD build kits`。

存圖（`kits_save_plot()`）另外需要：

* html：pandoc。RStudio 和 Quarto 都附了 pandoc；在 RStudio 外面用時，設 `Sys.setenv(RSTUDIO_PANDOC = "<pandoc 所在的資料夾>")`。
* png：webshot2，和一個 Chromium 系的瀏覽器（Chrome、Edge、Chromium）。

## 第一個例子

兩個世代、一條連結；可以直接貼上執行，最後是 `Checks: OK`。

```r
library(kits)
set.seed(1)
betas <- list(i1 = -2, i2 = -1, i3 = c(-0.5, 0.5), i4 = 0, i5 = c(0.5, 1.5), i6 = 2, i7 = 3)
young <- kits_sim(rnorm(500, 0), betas[1:5], prefix = "y")   # 0 資料：年幼世代答 i1–i5
old   <- kits_sim(rnorm(500, 2), betas[3:7], prefix = "o")   #   年長世代答 i3–i7；i3–i5 是共同題
item_y <- kits_item(young, scale = "young")                  # 1 校準：每個世代各一次
item_o <- kits_item(old, scale = "old")
link <- kits_linking(item_o, item_y)                         # 2 世代連結：old 的量尺 → young 的量尺
s_old <- kits_subj(old, item_o) |>                           # 4 計分與換算
  kits_convert(link)                                         #   old 先換到 young
s <- rbind(kits_subj(young, item_y), s_old)                  #   再疊
attr(s, "scale")                                             # "young"
kits_check(item_y, item_o, link)                             # 6 檢查
plot(link)
```

註解的數字是流程的步驟（見「流程」）；這個例子沒有 3 接到 D 和 5 常模。

怎麼讀：

* `item_y`：每一列是一題的一步；`thr` 是門檻。
* `link`：`a`、`b` 是一條直線；`dropped` 是篩掉的點。
* `s`：`theta` 是能力，`se` 是標準誤，`flag` 是 0 分、滿分或不適配。

下一步：

* 接到 D、做常模：入門篇第 3、5 節，`vignette("kits-flow", package = "kits")`。
* 加新的波次：進階篇第 8 節，`vignette("kits-advanced", package = "kits")`。
* 寫論文的方法段：`kits_algorithm(link)` 列出步驟、公式和文獻；可以複製的 LaTeX 在 `tutorial_algorithms/03_linking.html`。

## 速查表

| 步驟 | 我想要…… | 函數 | 最少要給 | 得到 | 注意 |
|---|---|---|---|---|---|
| 0 | 模擬作答 | `kits_sim(theta, betas, prefix)` | 能力、每題的步驟參數 | 矩陣 | 兩次呼叫的 id 不會撞名；要好讀的 id 時給 `prefix` |
| 0 | 作答表（選用） | `kits_resp(long)` | 長格式的作答＋孩子的欄 | 一列一個人次 | 一次一個領域；之後 `kits_item(resp)`、`kits_subj(resp, item)` 不用再對齊 |
| 1 | 一個世代一步做完 | `kits(X, scale = )` | 作答、量尺名 | `$item`、`$subj` | 不能直接 `kits_convert()`，用 `$subj` |
| 1 | 校準 | `kits_item(X, scale = )` | 作答、量尺名 | 門檻表（`coef()`） | 分數要從 0 開始、不跳號 |
| 2 | 世代連結 | `kits_linking(from, to)` | 兩個校準 | 直線 intercept + slope × θ | `from` 是要被換的，`to` 是目標 |
| 2 | 連結誤差 | `kits_linking(..., n_boot = 200)` | 同上 | `intercept_error`、`slope_error` | 每條都要有，才有 `se_link` |
| 2 | 連結一步一步 | `kits_link_points()` → `kits_link_screen()` → `kits_link_line()` → `kits_link_error()` | 同上 | 每一步的結果 | 除了 `info`（呼叫紀錄），和 `kits_linking()` 完全相同 |
| 3 | 接到 D | `kits_linking(item, kits_anchors())` | 校準、錨題表（`item_id`、`tau`） | 直線 | 至少 3 題錨題；一條連結只用一個 key |
| 4 | 能力 | `kits_subj(X, item, persons)` | 作答、校準 | 每人次一列 | `persons` 欄名：`visit_id`、`wave_m`、`age_m`、`sex`、`cohort` |
| 4 | 換量尺 | `x \|> kits_convert(link1) \|> kits_convert(link2)` | 分數或門檻、連結 | 同樣的表，新量尺 | 一行一條連結，照順序 |
| 4 | 疊表 | `rbind(a, b)` 或 `list_rbind(list(a, b))` | 同量尺的表 | 一張表 | 量尺不同會停；先換算，再疊；`bind_rows()` 會警告 |
| 4 | 原始分換算表 | `s \|> kits_subj_table()` | 分數（有 `wave_m`） | 每波次每原始分一列 | 只算答了完整題本的孩子；`filter()` 之後人數照樣對 |
| 4 | 每個作答組合的換算表 | `kits_subj_table(item)` | 校準 | 每個作答組合每原始分一列 | 不用資料 |
| 5 | 常模 | `kits_norm_screen()` → `kits_norms()` → `kits_norm_z()` | 分數＋`wave_m`、`age_m` | 百分位表；分數表多 `z`、`pct` | 不給 `ages` 時用 `kits_norm_ages()` 的規則；波次少時用 `degrees = c(3, 2, 1, 1)` |
| 5 | 照組的常模 | `kits_norms(s, by = "sex")` | 分數＋組的欄 | 一張表，組的欄在前 | 每組和只用這組的孩子配的常模相同；`plot(norm, scores = s)` 一組一格 |
| 5 | 常模一步一步 | `kits_norm_data()` → `kits_norm_fit()` → `kits_norm_table()` | 分數、月齡 | 每一步的結果 | 數字和 `kits_norms()` 相同；但沒有量尺標籤，`kits_norm_z()` 無法檢查量尺 |
| 6 | 題目適配 | `kits_itemfit(X, item)` | 作答、校準 | infit、outfit、flag | 預設是條件的（給定題本內的原始分），題數少時也不偏低；`type = "wle"` 在題數少時會小於 1 |
| 6 | 全部檢查 | `kits_check(...)` | 多個結果 | 警告和註記 | 疊過的表也可以給 |
| 6 | 來源 | `kits_info(x)` | 任何結果 | 呼叫、版本、md5 | |
| 6 | 一列的摘要 | `as.data.frame(x)` | 連結、校準、drift、info 等 | 一列（或每波次一列）的 data.frame | 幾個結果用 `bind_rows()` 疊起來比較 |
| 6 | 步驟在各組相同嗎 | `kits_dif(X, "raw")`、`kits_dif(resp, "wave")` | 作答、分組 | `print(d)` 印出 LR 檢定，`attr(d, "diff")` 是每步的 Wald z | `"raw"` 是 Andersen 的中位數分組；`plot(d, type = "steps")` |
| 7 | 性別 DIF | `kits_dif(X, sex, ref = )` | 作答、性別 | 每題的差 | 性別不是 boy/girl 時給 `ref`；作答表用 `kits_dif(resp, by = "sex")` |
| 7 | 題目表、成長曲線、階梯 | `kits_item_table(thr, norm, data = resp, dif = , stems = )` | D 上的門檻、常模 | 一題一列；給 `data` 時有通過月齡 `age50` | `plot(tab)` 畫成長曲線，`plot(tab, type = "ladder")` 畫階梯 |
| 8 | 偏移檢查 | `kits_drift(X_new, item_old)` | 新作答、舊校準（或只有 `item_id`、`item_score`、`beta` 的表） | `$fix`、`$items` | 新題不判斷；最高分比舊校準高的題不判斷，但建議固定 |
| 8 | 新波次固定參數 | `kits_item(X_new, fixed_params = drift$fix)` | `kits_drift()` 的結果 | 在原量尺上的校準 | 對**原本的**校準做 drift；一波接一波的題庫見 `?kits_item` |
| 9 | 演算法 | `kits_algorithm(x)` | 任何結果，或主題（`"cml"`、`"link"` 等） | 大約 30 行：步驟和公式、函數、文獻、教學和腳本 | 可以複製的 LaTeX 在 `tutorial_algorithms/` |
| 9 | 拆開的演算法 | `kits_stats(X) \|> algo_cml() \|> algo_thresholds()` | 數字 | 數字 | `algo_` 函數不記量尺、不檢查輸入 |
| — | 畫圖 | `plot(x)`、`plot(x, type = )`、`plot(x, items = )` | 校準、連結、常模、DIF、適配、題目表、偏移、篩波次等 | plotly 圖 | `type` 換一種圖（第一個是預設），`items` 選題號 |
| — | 存圖 | `kits_save_plot(p, "a.html", "a.png")` | plotly 圖 | 檔案 | html 要 pandoc，png 要 webshot2 |
| — | 信賴區間等 | `confint()`、`vcov()`、`logLik()`、`AIC()`、`nobs()`、`update()`、`simulate()` | 校準、連結、DIF、常模 | R 內建的結果 | `?kits_methods` |

## 流程

```text
0 資料    1 校準                       2 世代連結、3 接到 D
X03 ──kits_item()──▶ item03 ──kits_linking(item03, anchors)──▶ d_link：lan M03 → D
X36 ──kits_item()──▶ item36 ──kits_linking(item36, item03)───▶ co：    lan M36 → lan M03

4 計分與換算
kits_subj(X36, item36, persons) ──kits_convert(co)──▶ ──kits_convert(d_link)──▶ s36（D）

5 常模
kits_norm_screen() ──▶ kits_norms() ──▶ kits_norm_z()

6 檢查
kits_check(item03, item36, co, d_link, s36, ...)
```

M36 的分數要經過兩條連結，因為只有 M03 的量尺直接接到 D。

| 步驟 | 函數 | 在哪裡讀 |
|---|---|---|
| 0 資料 | 作答矩陣或 data.frame、`persons`；`kits_sim()`；`kits_resp()` | 入門篇第 0 節 |
| 1 校準 | `kits_item()`；一個世代時 `kits()` | 入門篇第 1 節 |
| 2 世代連結 | `kits_linking(item36, item03)` | 入門篇第 2 節 |
| 3 接到 D | `kits_linking(item03, anchors)`（`to` 是錨題表；KIT 的錨題 `kits_anchors()`） | 入門篇第 3 節 |
| 4 計分與換算 | `kits_subj()`、`kits_convert()`、`kits_subj_table()` | 入門篇第 4 節 |
| 5 常模 | `kits_norm_screen()`、`kits_norm_ages()`、`kits_norms()`、`kits_norm_z()` | 入門篇第 5 節 |
| 6 檢查 | `kits_check()`、`summary()`、`kits_itemfit()`、`kits_info()`、`plot()` | 入門篇第 6 節 |
| 7 題目分析與報告的圖 | `kits_dif()`、`kits_dif_combine()`、`kits_item_table(data = )` 和它的 `plot()` | 進階篇第 7 節 |
| 8 新的波次 | `kits_drift()`、`kits_item(fixed_params = )`、`plot(kits_vertical(s))` | 進階篇第 8 節 |
| 9 演算法和程式 | `kits_algorithm()`、`algo_` 函數、`system.file("scripts", package = "kits")` | 進階篇第 9 節 |
| 10 量尺的性質（Kolen & Brennan） | `kits_convert()` 的 `se_link`、`kits_vertical()`、`kits_precision()`、`kits_expected_score()`、`kits_link_check()`、`summary(link, theta = )` | 進階篇第 10 節 |

入門篇是 `vignette("kits-flow", package = "kits")`，進階篇是 `vignette("kits-advanced", package = "kits")`。
每個函數的說明頁（例如 `?kits_linking`）都有可以直接跑的例子。

## 名詞表

| 中文 | English | 欄名／函數 | 一句話 |
|---|---|---|---|
| 部分給分模型 | PCM (partial credit model) | `kits_item()` | 每題的分數是 0、1、…、m，每一步有一個步驟參數；只有 0、1 分時就是 Rasch 模型。 |
| 步驟參數 | step parameter | `beta` | PCM 裡從 k − 1 分到 k 分的參數；同一題的 `beta` 不一定由小到大。 |
| 條件最大概似 | CML (conditional maximum likelihood) | `kits_item()` | 以原始分為條件估題目參數，不用假設孩子能力的分佈。 |
| ENORM | ENORM | — | dexter 對同一個模型的名字。 |
| 能力、logit | ability, logit | `theta` | 孩子在量尺上的位置；一個校準自己的量尺用 logit 當單位。 |
| 加權概似估計 | WLE (weighted likelihood estimate) | `kits_subj()` 的 `theta` | 每個孩子的能力估計；0 分和滿分也有有限的值。對 PCM，它就是 Jeffreys 眾數估計：log L + ½ log I 的最大值。 |
| WLE 的標準誤 | SE of the WLE | `se` | 1/√I(θ̂)，I 是測驗訊息（Warm, 1989 的漸近 SE）。 |
| 門檻 | threshold | `thr`（公式寫成 thr_ik） | 能力在這裡時，得到 k 分以上的機率是 .5。 |
| Thurstonian 門檻 | Thurstonian threshold | `thr` | 門檻在文獻上的名字：累積機率（k 分以上）是 .5 的能力。 |
| RP50 | response probability 50 | `coef(item, rp = )` | 門檻用的機率是 .5；RP80 用 .8，題目會放得更高。 |
| milestone step | milestone step | `item_score == 1`、`step = 1` | 每題的第 1 步，也就是「通過」（分數至少 1）；通過一題就像達到一個發展里程碑。 |
| 世代 | cohort: M03, M36 | `cohort` | KIT 的兩個樣本：從 3 個月開始追蹤，和從 36 個月開始追蹤。 |
| 波次、月齡 | wave, age in months | `wave_m`、`age_m` | 波次是預定的施測月齡，月齡是實際的月齡；`_m` 是月。 |
| 人次 | person-wave | 表的一列 | 一個孩子在一個波次是一個人次。 |
| 作答組合 | item pattern | `pattern_id` | 一個孩子實際答的那組題（哪幾題有作答）。 |
| 題本 | standard form（簡稱 form） | `kits_subj_table()`、`form_pct` | 一個波次大多數孩子答的那組題。 |
| 量尺標籤 | scale label | `scale`、`attr(x, "scale")` | 每個校準定義一條量尺，標籤跟著結果走，換算時用來檢查。 |
| 世代連結 | cohort link | `co`、`kits_linking(item36, item03)` | 用共同題把一個世代的量尺換到另一個世代；變數習慣叫 `co`。 |
| mean/sigma | mean/sigma | `kits_linking()` | 用共同題門檻的平均和 SD 求一條直線 intercept + slope × θ。 |
| 錨題 | anchor item | `anchors`（`item_id`、`tau`）；`kits_anchors()` | 和 GSED 相同的題，用來接到 D。 |
| GSED | Global Scales for Early Development | — | WHO 的嬰幼兒發展量表，D 分數的來源。 |
| D 分數 | D-score | 換算後的 `theta` | GSED 的發展分數；接到 D，KIT 的分數就可以和其他國家比。 |
| GSED 難度 | tau | 錨題表的 `tau` | 錨題在 GSED 的難度，單位是 D；它不是 KIT 的門檻 `thr`。 |
| shash | sinh-arcsinh distribution | `kits_norms()` 的 `mu`、`sigma`、`epsilon`、`delta` | 可以偏態、厚尾的分佈；四個參數各是月齡的多項式。 |
| 同齡 z | same-age z | `kits_norm_z()` | 孩子在同月齡常模裡的位置；百分等級是 `100 * pnorm(z)`。 |
| md5 | md5 | `kits_info()` | 資料的指紋：資料改了，md5 就不同。 |
| Lord-Wingersky 遞迴 | Lord-Wingersky recursion | `kits_precision()`、`algo_lw()` | 給定能力時，算原始分分佈的方法。 |
| 演算法函數 | algorithm functions | `algo_cml()` 等 12 個 | 演算法本身：給數字、回數字，不記量尺；`kits_` 的函數呼叫它們。 |

## 常見錯誤

* **分數從 1 開始**（例如 1、2、3）：kits 要分數從 0 開始。先減 1：`X <- X - 1`。
* **分數中間有缺**（例如只有 0 和 2）：先重新編碼成 0、1。
* **錨題表的欄名**：錨題表要 `item_id` 和 `tau`；門檻表才是 `item_id`、`item_score`、`thr`。
* **連結的方向**：`kits_linking(from, to)` 的 `from` 是要被換算的校準，`to` 是目標；給反了不會停，看印出來的 `scales:` 那一行。
* **量尺不同**（`kits_error_scale`）：先用 `kits_convert()` 換到同一個量尺，再疊表；連結照順序給。
* **`persons` 的欄名**：要叫 `visit_id`、`wave_m`、`age_m`、`sex`、`cohort`；別的名字會被忽略，到常模那一步才停下來。
* **沒給 `scale`**：量尺用資料物件的名字；資料不是一個有名字的物件時，叫 `"scale 1"`，並說一聲。一定要給 `scale`。
* **性別不是 boy、girl**：`kits_dif(X, sex, ref = "M")`；錯誤訊息會列出資料裡的值。
* **常模不收斂、百分位亂跳**：波次少時降低多項式的次數，`kits_norms(..., degrees = c(3, 2, 1, 1))`。
  `nlminb()` 差一點收斂就停下時，常模和 `kits_item()` 一樣，在自由的係數上再走幾步 Newton；已經收斂的結果不變。
* **存 html 失敗**：沒有 pandoc；設 `RSTUDIO_PANDOC`（見「安裝」）。
* **訊息沒有列出 id**：訊息和錯誤只給數目，不給 person id；要看 id 時，設 `options(kits.show_ids = TRUE)`。
* **`Checks: 1 note`**：不一定是錯；`summary(x)` 看是哪一項。`note` 是知道就好，`warn` 要處理或在報告中說明。
* **`nlminb code 1 (singular convergence)`**：狀態是 `ok` 時是正常的；它的停止標準很嚴，梯度已經夠小。
* **某一波的信度很低或是負的**：通常是題本對這一波太容易或太難；看 `summary()` 裡 `perfect` 和 `zero` 的比例。

## 我的資料不是 KIT

* **分數從 1 開始**：先減 1。
* **性別的寫法不同**：`kits_dif(..., ref = )` 指定參照組。
  `kits_norms(s, ages, by = "sex")` 照 `sex` 欄的值分組，代碼用什麼都可以。
* **世代不叫 M03、M36**：`kits_dif_combine(d1, d2, cohorts = c("A", "B"))`（沒給時用參數的名字）。
  `plot(v, type = "median", cohorts = c("A", "B"))` 選要畫的世代；`v` 是 `kits_vertical()` 的結果。
* **只有一個世代**：不需要連結；`kits(X, scale = "my scale")` 一步做完校準和計分。
* **沒有錨題**：跳過「接到 D」，留在 logit 的量尺上；常模的函數也可以用 logit 的分數（會說一聲分數不在 D 上）。
* **波次少**：`kits_norms(..., degrees = c(3, 2, 1, 1))`；預設的 `c(5, 4, 3, 2)` 是給 KIT 的 3 到 72 個月、十幾個波次。
* **月齡範圍不同**：`kits_norm_ages()` 從資料定範圍，不用改；圖的刻度用 `age_ticks`（預設是匯出的常數 `AGE_TICKS`）。
* **要別的百分位**：`kits_norms(..., pct = )`、`kits_norm_table(..., pct = )`；預設是匯出的常數 `NORM_PCT`（3、5、10、…、97）。

## 細節

### 資料：直接給

像 mirt、TAM 一樣，資料直接交給函數；模型是固定的（PCM），不用先包裝。
每個收作答的函數都收三種形式：

* 矩陣：一列一個人次，列名是 visit_id；一欄一題；NA 是沒施測。
* 寬的 data.frame：同樣的形狀；visit_id 是叫 `visit_id` 的欄（沒有時用列名），其他欄都是題。
* 長的 data.frame：`visit_id`、`item_id`、`item_score`。

人次、題號重複，分數不是 0 以上的整數，分數中間有缺，都會停下來（`kits_error_data`）。
作答組合分成幾群、群和群之間沒有共同題時，也會停下來：這樣的資料不在同一把尺上。
沒有作答的題和人次拿掉，並說一聲。

沒有 visit_id（矩陣沒有列名、data.frame 沒有 `visit_id` 欄）時，照列的順序編號 1、2、3……，並說一聲。

編號的欄（0.23.0 改名）：

| 欄 | 意思 | 0.22.0 以前 |
|---|---|---|
| `visit_id` | 一個孩子的一次施測（人次）；長格式裡同一個人次有很多列 | `person_id` |
| `subj_id` | 孩子；每一波都一樣 | `pid` |
| `pattern_id` | 作答組合（item pattern）：一個人次答了哪幾題 | `booklet_id` |
| `raw`、`raw_max` | 原始分、這個作答組合的滿分 | `booklet_score`、`booklet_max` |

資料裡還是舊的欄名 `person_id`、`pid` 時，這一版照常認得，並說一聲（`kits_message_data`）；下一版拿掉。

人的資料用參數給：

* `kits_subj(X, item, persons = persons)`：`persons` 有 `visit_id`，和 `subj_id`、`wave_m`、`age_m`、`sex`、`cohort` 中的幾欄，照 `visit_id` 接上。
  `subj_id` 是孩子的編號，每一波都一樣。
  其他的欄不留；要用時自己 `left_join()`。
  之後常模的函數直接用：`kits_norm_screen(s)`、`kits_norms(s)`、`kits_norm_z(s, norm)`。
* `kits_item_table(thr, norm, data = )`：`data` 是有 `age_m`、`wave_m` 的作答；分數至少 1 算通過，`age50` 欄是 50% 通過的月齡。
  每題的通過不同時用 `step =`，例如反向的 0–4 題至少 3 分算通過：`step = c(cogD01 = 3)`（沒列的題用 1）。
  每題實際的通過率對 θ̂（一組一條線，灰線是模型）：`plot(item, type = "observed", data = X, group = "wave")`。
* `kits_dif(X, sex)`：每一列一個性別。

也可以先把作答和孩子的資料放進一張表（選用）：

```r
resp <- kitm03$responses |>
  dplyr::filter(domain == "lan") |>
  kits_resp()
item <- kits_item(resp)              # 量尺叫 "lan M03"（領域和世代）
s <- kits_subj(resp, item)           # 孩子的欄跟著來
```

* `kits_resp()` 收長格式（一列一個作答，加上孩子的欄），或寬格式加上 `persons =`。
* 一張表只放一個領域。
* `kits_item_table(thr, norm, data = resp)`、`kits_dif(resp, by = "sex")`、`kits_itemfit(resp, item)` 也直接收它。
* `print(resp)` 有各波的人數，`summary(resp)` 有各波的題本，`plot(resp)` 畫每一波答了哪些題。

### 連結：`kits_linking()` 的兩種用法

| 用法 | `to` 是 | 用哪些點 | 說明頁的名字 |
|---|---|---|---|
| 世代連結 | 一個校準（`kits_item()`） | 共同題的所有步驟；`step = 1` 只用 milestone step | between two calibrations |
| 接到 D | 錨題表（`item_id`、`tau`） | 每題錨題的通過步驟（預設 `step = 1`） | to anchors |

連結誤差預設不算；要看時給 `n_boot`（例如 200；報告用 2000）。
重抽的種子預設是 `seed = 1`，所以每次結果都一樣。
同一個世代裡比孩子時，大家用同一條線，不需要連結誤差；跨世代比較或報告 D 的不確定性時才需要。
`from` 可以是有名字的 list，一個領域一個校準；這時幾個領域共用一條線。
沒有錨題的領域（例如社會情緒）也用這條線；印出來時會寫 `same line, no anchors: social`，`kits_info()` 也記下來。
這些領域的 KIT 量尺分數是 D 的單位，不是 GSED 的 D。

KIT 的錨題表在 `kits_anchors()`：29 題，`tau` 來自 dscore 的 key gsed2406（預設），`note` 記下 KIT 怎麼對題。
`domain` 欄是錨題的領域（lan、motor、cog），從 GSED 的題號來。
`kits_anchors("gsed2510")` 只有 8 題（2 題語言、6 題動作、沒有認知），D 的尺度也不同，所以會警告；兩個 key 的 D 不能互換。
一條連結只用一個 key：兩個 key 疊在一起時，`kits_linking()` 會停下來；連結會記下用的 key。

### 量尺：`scale` 和 `kits_convert()`

每個校準定義一條量尺，標籤是 `scale`（例如 `"lan M03"`）。
標籤會跟著結果走：`kits_subj()`、門檻表、`drift$fix` 都記得自己的量尺。
連結記得兩端：`link$from_scale`、`link$to_scale`。
`kits_convert()` 換量尺，用管線一行一條連結，照套用的順序寫：

```r
s36 <- kits_subj(X36, item36, persons = persons) |>
  kits_convert(co) |>        # lan M36 -> lan M03
  kits_convert(d_link)       # lan M03 -> D
```

`kits_convert(s36, co, d_link)` 是同一件事寫成一次呼叫，結果相同。
連結接錯、順序錯，會停下來（`kits_error_scale`）。
每條連結都是一條直線：新的值 = intercept + slope × 舊的值（SE 乘上 slope）。
intercept 是截距，slope 是斜率。

**kits 0.23.0 起不用 a、b 這兩個名字。**
同樣的字母在不同的地方意思不同：

| 來源 | 截距 | 斜率 |
|---|---|---|
| kits 0.22.0 以前、kit_five 的連結表 `02_link_table.csv` | `a` | `b` |
| Kolen & Brennan、TAM、plink、equateIRT | B（或 b） | A（或 a） |
| kits 0.23.0 | `intercept` | `slope` |

其他的名字：`a_error`、`b_error` 改叫 `intercept_error`、`slope_error`；`a0`、`b0` 改叫 `intercept0`、`slope0`；
`draws` 的欄叫 `intercept`、`slope`；`ab_cov` 是 `vcov`（2 × 2 的矩陣）的對角線外的那一格。
數字都不變。
舊的程式寫 `link$a`、`link$b` 時會停下來，說新的名字（0.24.0 拿掉這個檢查）。
`coef(link)` 給 `c(intercept = , slope = )`。

### 疊表和 dplyr

結果的表都是 data.frame（不是 tibble），dplyr 的動詞直接可以用。
`filter()`、`mutate()`、`arrange()` 之後，類別和量尺留著，描述整張表的診斷和來源拿掉（和 `x[i, ]` 相同）。
疊兩個 kits 的表用 `rbind()` 或 `purrr::list_rbind()`。
兩個都會檢查量尺，也會把兩張表的紀錄（換算、原始分對照）合在一起。
量尺不同時會停下來（`kits_error_scale`），所以先換算，再疊起來。
`dplyr::bind_rows()` 只留第一張表的紀錄，所以疊有紀錄的 kits 表時會警告（`kits_warning_scale`）。
和普通的 data.frame 疊在一起時，會警告（`kits_warning_scale`）。
疊好的表，`summary()` 和 `kits_check()` 照樣能用：`kits_subj()` 的表會重算診斷，其他的表給一個註記。

### 印出來的表

表格印小數點下 2 位；物件裡的數字沒有被捨入。
超過 30 列的表只印前 10 列，最後一行說還有幾列。
`print(x, n = 50)` 印更多列（校準也可以）；`as.data.frame(x)` 拿到全部的列。
固定參數的題，印出來會標出來（SE 是 0）。

### 檢查：`summary()` 和 `kits_check()`

每個結果都會檢查自己。
印出來的最後一行是 `Checks: OK`，或有幾個警告、幾個註記。
`summary(x)` 列出所有檢查。
`kits_check(...)` 把好幾個結果的警告和註記放在同一張表；沒有診斷的物件也列成一列，不會停下來。

### 來源和隱私：`kits_info()`

每個結果記下呼叫、kits 版本、設定、資料的指紋（md5）。
呼叫只記下名字；資料直接當值給時（例如用 `do.call()`），記成 `<data>`。
沒有時間戳記，所以同樣的輸入給完全相同的物件。
訊息和錯誤只給數目，不列出 person id；`options(kits.show_ids = TRUE)` 會列出來。

### 演算法和程式：`kits_algorithm()`、`algo_` 函數

`kits_algorithm(x)` 寫出結果背後的演算法，大約 30 行：
步驟，一步一個公式；跑這些步驟的 `kits_` 和 `algo_` 函數，用 `|>` 接起來；2 到 4 篇文獻；去哪裡讀更多。
給結果時，`This result` 一行是這個結果自己的數字和量尺，例如連結的線和「from -> to」。
也可以給主題：`"cml"`、`"wle"`、`"link"`、`"norms"`、`"norm_screen"`、`"dif"`、`"drift"`、`"information"`。

`algo_` 是 12 個演算法函數：給數字、回數字，不記量尺、不檢查輸入；`kits_` 的函數呼叫它們。

* 校準：`algo_pcm()`、`algo_esf()`、`algo_cml()`、`algo_thresholds()`；分數：`algo_jme()`。
* 連結：`algo_screen()`、`algo_mean_sigma()`、`algo_link_error()`、`algo_line()`。
* 常模：`algo_shash()`、`algo_shash_z()`；原始分的分佈：`algo_lw()`。

分數還是叫 WLE。
對 PCM，WLE 就是 Jeffreys 眾數估計 argmax log L + ½ log I（Warm, 1989; Firth, 1993; Magis, 2015），kits 用 RTMB 算。
有兩個高點時取比較高的那一個；SE 是 1/√I(θ̂)。

`kits_algorithm()` 的最後一行指到兩個地方：

* 教學（`tutorial_algorithms/` 的六份：`01_calibration`、`02_scoring`、`03_linking`、`04_norms`、`05_checks`、`06_scale`）：
  每個公式都有排好的樣子和可以複製的 LaTeX，寫論文的方法段用；管線貼到 R 就能重算。
* 教學腳本和它的段落（`system.file("scripts", package = "kits")`：`cml.R`、`wle.R`、`link.R`、`norm_screen.R`、`norms.R`、`dif.R`、`drift.R`、`information.R`）。
  它在示範資料上不用 kits 算一次（RTMB、base R、dplyr、tidyr、purrr），最後和 kits 的結果比。
  第一段 `inputs` 換成自己的資料就能用。

```r
kits_algorithm(item_y)                     # 校準：步驟、公式、文獻
kits_algorithm(link)                       # 這條線和「from -> to」
kits_algorithm(s)                          # 分數：WLE 和個人適配
kits_algorithm("norms")                    # 給主題
kits_stats(young) |> algo_cml() |> algo_thresholds()   # item_y 的門檻
system.file("scripts", "cml.R", package = "kits")
```

校準也包著兩個 θ 的函數：`kits_information(item_y)`（測驗訊息）和 `kits_expected_score(item_y)`（期望原始分）。
`kits_pcm_prob(item_y, theta)` 把它們一題一題拆開：類別機率 P(X = x | θ)、期望分數、變異數（題目訊息）；`kits_esf()` 是給定原始分 r 的另一面。

`kits_inter()`（交互作用模型）0.23.1 拿掉了：題目適配看 `kits_item()` 的條件 infit、outfit（`item$itemfit`、`kits_itemfit()`）；
dexter 的交互作用模型在 `kitsPreview::kits_inter()`。

### 拿掉的函數和寫法（0.23.0、0.25.0）

0.22.0 停用的 11 個函數和幾個寫法，0.23.0 拿掉了；0.25.0 拿掉的在表的最後三列。
0.18.0 以前的參數名字 `dataSrc`、`D` 也拿掉了：現在叫 `data`、`d`。

| 拿掉的寫法 | 現在的寫法 |
|---|---|
| `kits_age50(X, age_m, wave_m)` | `kits_item_table(thr, norm, data = )` 的 `age50` 欄 |
| `kits_mean_sigma(from, to, id)` | `points \|> kits_link_screen() \|> kits_link_line() \|> kits_link_error()` |
| `kits_link_se(co, d_link, theta = x)` | `tibble(theta = x) \|> kits_convert(co, d_link, from = "lan M36")` |
| `kits_link_screening(link, theta)` | `summary(link, theta = )$screening` |
| `kits_wle_table(item)` | `kits_subj_table(item)` |
| `kits_subj_table(item, X, wave_m)` | `kits_subj(X, item, persons) \|> kits_subj_table()` |
| `kits_plot_cohort(自己做的表, scores)` | `plot(kits_vertical(s), type = "median")` |
| `kits_plot_norms(自己疊的表, scores)` | `kits_norms(s, by = "sex") \|> plot(scores = s)` |
| `kits_item_positions()`、`kits_ladder_marks()` | `kits_item_table(thr, norm, age_pass = , dif = , stems = )` |
| `kits_plot_growth(norm, items)` | `plot(kits_item_table(...))` |
| `kits_plot_ladder(items, marks)` | `plot(kits_item_table(...), type = "ladder")` |
| `s$z <- kits_norm_z(s, norm)` | `s <- s \|> kits_norm_z(norm)`（多 `z`、`pct` 兩欄；舊的寫法會停下來說清楚） |
| `plot(norm, d = s$theta, age_m = s$age_m)` | `plot(norm, scores = s)` |
| `kits_item(dataSrc = X)`、`kits_norms(D = d, ...)` | `kits_item(data = X)`、`kits_norms(d = d, ...)` |
| `kits_code(x)`、`kits_equations(x)`、`kits_rtmb_code(x)` | `kits_algorithm(x)`；可以複製的 LaTeX 在 `tutorial_algorithms/` |
| `kits_lrtest(X, item, split = "raw")` | `kits_dif(X, "raw")`（印出 LR 檢定）（也可以照 `"wave"`、性別分組） |
| `kits_age_pass(X, age_m, wave_m)` | `kits_item_table(thr, norm, data = )` 的 `age50` 欄 |

`points` 是自己準備的點：一個 data.frame，有 `id`、`item`、`from`、`to` 四欄。
點的量尺用 `kits_link_screen(points, from_scale = , to_scale = )` 給。

`kits_norms(s, ages)` 和 `kits_norm_z(s, norm)` 照第一個參數的類別（S3）分：給 `kits_subj()` 的結果時，月齡來自它的 `age_m` 欄。
教學腳本 `norms.R` 裡，分數的變數照公式叫 `D`。

### 畫圖（0.22.0）

每個結果用 `plot()` 畫。
`type =` 換一種圖，第一個是預設。
`items =` 選要畫的題，給題號。
新的圖有：校準的類別機率曲線和 Wright 圖、性別 DIF、題目適配、偏移、篩波次、作答表、題目表。
照組的常模一組畫一格；給 `scores =` 時，每一格畫自己那一組的孩子。
多的資料用名字給，例如 `plot(norm, scores = s)`、`plot(item, type = "wright", scores = s)`。

### R 內建的方法（0.22.0）

R 使用者已經會的名字，不用記新的函數（`?kits_methods`）：

* `confint()`：校準的門檻、連結的 intercept 和 slope（有重抽時用重抽的分位數）、DIF。
* `coef(link)`：連結的 `c(intercept = , slope = )`，順序和 `lm()` 一樣。
* `vcov()`：校準的步驟參數的共變異數矩陣。
* `logLik()`、`nobs()`，所以也有 `AIC()`、`BIC()`：校準和常模；只比較同一份資料的配法。
* `update()`：改一個設定重配，例如 `update(norm, degrees = c(4, 3, 2, 1))`。
* `simulate(item, theta = )`：從校準的參數模擬作答。

步驟函數的說明合成三頁：`?kits_item_steps`（校準）、`?kits_link_steps`（連結）、`?kits_norm_steps`（常模）。

### 錯誤、警告、訊息

每個錯誤、警告、訊息都有類別，可以用類別接住，例如 `tryCatch(..., kits_error_scale = function(e) ...)`。
類別的清單在 `?"kits-package"`。

## 進階

### 性別 DIF：`kits_dif()`、`kits_dif_combine()`

方法和 dexter 的 `DIF()` 相同（Bechger & Maris, 2015），但不用 dexter：男女各自用 kits 的 CML 校準，再比較相對的題目參數。
`ref = "boy"` 時，正的有利女孩。
`group` 也可以是 `"raw"`（Andersen 的中位數分組）、`"wave"` 或欄名：`print(d)` 印出步驟在各組相同的 LR 檢定，`attr(d, "diff")` 是每一步的 Wald z，`plot(d, type = "steps")` 畫出來。
`kits_dif_combine(d1, d2)` 合併兩個世代；世代的名字用 `cohorts` 給，沒給時用參數的名字。
`DIF` 欄給成長曲線和階梯的圖上色：沒有灰色、有利女孩橘色、有利男孩藍色、不一致黑色。

### 量尺的性質（Kolen & Brennan, 2014）

這一批方法來自 Kolen 和 Brennan（2014）的 *Test Equating, Scaling, and Linking*（第 3 版）。
例子和頁碼在進階篇第 10 節。

* **每個 θ 的連結誤差**：`n_boot > 0` 的連結留下重抽的線。
  每條連結都有重抽時，`kits_convert()` 加上 `se_link` 和 `se_total`；有一條沒有時，不加欄。
  選定的 θ 也這樣算：`tibble(theta = x) |> kits_convert(co, d_link, from = "lan M36")`。
  一串連結的重抽是獨立的，用動差直接算。
* **斜率和 1 比**：PCM 的連結只需要平移；kits 保留 mean/sigma 的自由斜率，斜率和 1 差太多時有註記。
* **接到 D 的用語**：用錨題接到 D 是校準（量尺對齊），不是等化；寫方法段時照這個用語（進階篇第 10.3 節）。
* **垂直量尺**：`kits_vertical(s)` 每個世代 × 波次一列：平均、SD、中位數、P25、P75、扣掉測量誤差的 `sd_true`、成長、效果量、`p_sup`。
  `plot(v)` 畫平均 ± SD，`plot(v, type = "median")` 畫中位數和 P25 到 P75。
  好幾個領域時，用 `list_rbind(names_to = "domain")` 疊起來，一個領域畫一格。
  0 分或滿分很多的波次，`sd_true` 要當成下限來讀。
* **精準度**：`kits_precision(item, forms)` 看每個題本在量尺的每個位置量得多準（偏誤、SD、RMSE）；原始分的分佈用 Lord-Wingersky 遞迴。
* **原始分對照**：`kits_subj_table(item)` 的 θ 代進另一組題的 `kits_expected_score(item, items = )`；PCM 的原始分是充分統計量，原始分 → θ 的表就是計分（§6.4）。
* **世代連結的獨立檢查**：`kits_link_check()` 用同齡的兩群和合併校準（concurrent CML），各估一次同一條線，再比較。
  同齡線的收斂和 `kits_item()` 一樣看梯度和 Hessian。
  兩個世代的 person id 重複時，用 `persons_from`、`persons_to`，一個世代一張表。
* **篩點的規則**：`kits_linking()` 的 `spread`（`"sd"`、`"mad"`、`"se"`）和 `unit`（`"step"`、`"item"`）。
  `summary(link, theta = )` 看篩點把換算值移動多少。
  `"se"` 只適合小樣本：KIT 的 `thr_se` 很小，篩掉的點會太多。
* **反應機率**：`coef(item, rp = 0.8)` 給 RP80 的門檻；`kits_item_table(..., rp = 0.8)` 給 80% 的孩子通過的月齡（`age80`）；接到 D 一直用 RP50。
  完全分離的題（小的孩子都沒通過、大的都通過）列在題目表的診斷裡，月齡只知道在 `age_last_fail` 和 `age_first_pass` 之間。
* **品質管制**：作答組合太小、共同題比例太低、錨題集中在一個領域、兩個世代同齡的結果不一致時，都有註記。
  新波次的檢查表在 `?kits_check` 和進階篇第 10.10 節。

`kits_vertical()`、`kits_link_check()` 沒有 `kits_algorithm()`；`kits_precision()` 的是 `"information"`。

### 和 kitsPreview 的數字

姊妹套件 kitsPreview 用 dexter 和 cNORM 算同樣的結果，用來驗證；新手不需要它。
kits 的 PCM 就是 dexter 的 ENORM（分數 0..m）。
同樣的模擬資料（套件的測試）：

* β、門檻、θ：差距在 1e-6 以內；β 和門檻的 SE 也是。
* θ 的 SE：kits 0.25.0 起用 1/√I，dexter 用 Warm 的式子，所以 kits 的稍小（中間的分數約小 0.2%，0 分和滿分約小 20%）。
* infit、outfit（kits 的 `type = "wle"`；kitsPreview 只有這一種）：1e-5 以內。
* 連結、篩波次、月齡範圍：完全相同。
* 常模的 P 欄：1e-4 以內。
* 性別 DIF：`dif` 在 1e-6 以內；`kits_dif_combine()` 的判斷完全相同。
  整體檢定在每題都留著時相同；有題目被拿掉時，kits 的檢定不受原點影響，和 dexter 的 `DIF()` 不同。
