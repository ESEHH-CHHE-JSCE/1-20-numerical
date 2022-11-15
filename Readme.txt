例題1.20サンプルコードについて

■添付一覧
　このサンプルには以下のファイルが含まれます．

・「common.h」一部の変数を指定したヘッダファイル
・「condition.d」計算条件設定ファイル
・「boundary.d」境界条件設定フラグ記述ファイル
・「grid.f」計算格子作成プログラムソースコード
・「density2d.f」鉛直2次元密度流解析プログラムソースコード
・「csvdata.f」


■使用の手順

１）上記の「grid.f」をコンパイルして実行することで，「common.h」，「condition.d」，および「boundary.d」に記述された条件に応じて計算格子が作成される．
計算条件に用いた変数については下記のとおりである．

◇「common.h」設定変数について
im,jm：x方向，z方向に用意する配列数，計算格子の数に境界条件用の配列をプラスした数

◇「condition.d」設定変数について
nx,nz：x方向およびz方向格子数
iscf：ソース箇所設定の有無
chleng,height：流路長および流路高さ
g,rho,rnu,diffyr：重力加速度，周辺流体の密度，動粘性係数，拡散係数/動粘性係数の比
tuk,etime,dt：計算結果出力時間刻み，計算終了時刻，計算時間刻み
lsor,soralp：圧力の計算の反復数，圧力計算の緩和係数
pretime：初期条件からの緩和にかかる時間

◇「boundary.d」設定変数について（原則として変更不要）
ib1u,ib1w,ib2u,ib2w：上流および下流境界の流速成分の境界条件設定用フラグ
ib3u,ib3w,ib4u,ib4w：上面および下面境界の流速成分の境界条件設定用フラグ
ib1p,ib2p,ib3p,ib4p：圧力計算の境界条件設定用フラグ
ib1r,ib2r,ib3r,ib4r：濃度計算の境界条件設定用フラグ

◇「grid.f」設定変数について
・84行目から96行目
c     source settings
      do i=nx/2,nx/2
        do k=2,4
            isw_w(i,k)=1
            scw_w(i,k)=0.7
        end do
      end do
      do i=nx/2,nx/2
        do k=3,5
            isr_c(i,k)=1
            scr_c(i,k)=rho+5.0
        end do
      end do

　密度噴流ソースの位置と流速を決定

２）その後，「density2d.f」をコンパイルして実行することで，計算結果を得る．

３）この時に生成されたファイル群に対し，「csvdata.f」をコンパイルして実行することで，出力時間刻みで生成された濃度および流速ベクトル成分のcsv形式の結果ファイル群が作成される．


■動作条件

これらコードはWindows10Pro上のIntel oneAPI 2021による動作を確認している．
