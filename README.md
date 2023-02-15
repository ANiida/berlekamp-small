# berlekamp-small

lyon.c：　バーレカンプマッシー復号

p-tar.c:　ピーターソン復号

# 20230215

バーレカンプマッシー法を使った、Niederreiter暗号が完成しました！

根気よく手直ししてくださったるばーとさんに特別な感謝を申し上げたい。

# 20230206

結局自分が必要としていたのはIPsecとかVPNだったんだということに気づく。

もうこんなに必死になって暗号をやる必要はないんだと安心する（のはまだ早いかも）。

もっと楽しんでプログラムできるようなテーマをやってみたい。
暗号の方もダブルサイズのシンドロームと、３の拡大体上のGoppa符号をやめればほぼ開発終了。
いろんな復号法がある中で、それぞれがどのような関係にあるのかとか、暗号目的でない研究ができそう。

インターネットを監視から救う方法は整っている。
IPｖ６に組み込まれている暗号化機能を設定すれば、誰でも必要な通信経路を暗号化できる。
そしてそのうち特定の拠点間だけでなく、全てのルーターが暗号化と改ざんを防げるようになれば、
私の夢見た暗号化ネットワークは実現するのだ。

みんな知らないだけw

# 20230205

やっと位置づけがはっきりしてきました。

berlekamp-smallはbm法で復号するリポジトリ。
GoppaDecorderはパターソンとユークリッド方の両方で復号するリポジトリ。

もうすぐ完成です。

色々確かめていくうちに、今までずっとリードソロモンで復号していたことが明らかになって、急遽修正に。
で、やっと目処が付きました。
確認は大事。

あと他の人が読んでもわかるようにRubatoさんの指示に従ってメンテナンスに移行します。

ピーターソンとバーレカンプマッシーは、シンドロームだけで復号できるので、行列の値とかは気にしなくてすむのに、
ユークリッドとパターソンはゴッパ多項式を必要とするので、なぜ同じものを復号するのにゴッパ多項式が必要なのだろうと疑問に思う。

情報が少ないほうが有利であるというのなら、速度的な面からもｂｍ法に軍配が上がる。

あと、
$GF(2^m)$
の２乗の線形性という性質を使って、半分のシンドロームから２倍サイズのシンドロームを生成する実験をする予定。

オックスフォード出版局の出している符号の本で、バイナリGoppaの性質というところでほんの少し触れているだけなのだが、
これでダブルサイズシンドロームの謎が解けるかもしれない。


# 20230123

なんとワイルド系符号も、２次元巡回符号も、シンドローム復号も全部日本人が重要な業績を残しているんですね。
これを合体させたら見たこともない暗号になるのかも？！
ワイルド２次元シンドローム復号問題！ｗ

# 20230122

ピーターソン復号法を実装。
２，３時間で終わると思ってたら５時間超えた。

狙いは３以上の拡大体で構成したリードソロモン符号を、ピーターソン復号すること。
代数幾何符号は難しいので奇標数の拡大体の方にナイーブに一般化する予定。
まだ狙ってるだけで不可能であることがわかるかもしれないけど、もし鍵方程式を解くやり方と関係していれば面白いことになる。

つまりワイルドリードソロモンのわかりやすいシンドローム復号法ができるということだ。
復号法さえできれば速度の違いこそあれ、ワイルドマックエリスを動かすことができる。

ユークリッドは微分しないといけないし、バーレカンプは動くかどうか怪しいし、パターソンは複雑で難しいし、同じシンドローム復号の中でも遅いけど一番ラクな方法がピーターソン。

シンドロームを並べてガウスの消去法で連立方程式を解くだけ！
ピーターソンはすごいなあ。

遅くてもいいからなるべく短時間で結果を出したかっただけ。
とてもわかり易い復号法。

リスト復号だとシンドローム復号問題を使うことができないのでとりあえず冒険しない方針で実装。

シンドローム復号と呼ばれるものは全て同じデータだけに依存している。
つまり各々の復号法は互いに関連しあっていて、違うように見えるだけで同じなのかもしれない。

ピーターソン凄すぎ。
Goppaでも復号できちゃう。
ケンブリッジの学部向けの符号の教科書で実装していたから、最初はGoppaを復号するためにはユークリッドでなければならないという思い込みがあった。
一体今までの数年間の苦しみは何だったのか？ｗ

5時間で解決できる問題を、大学を卒業できる年数だけ遠回りして到達したのだ。
今の所マックエリスは安全という前提で作っているけど、Goppaには代数的構造の情報が多く含まれているので、そのうちスクランブラでも隠しきれない攻撃が出てくる可能性がある。
次のターゲットとしては、符号の同一性判定問題とかだけど、色々符号化方法を集めるだけでなく、やはり安全安心の暗号を目指さなければならない。

Known Syndromeだけでフェンラオ限界まで誤り訂正できる復号法があれば素敵じゃない？
夢物語かもしれないけど・・・

ピーターソンは簡単だった。
あとは今井センセの符号理論でｋｗｓｋ動作原理を理解するだけ！


# 20230121

今、mine以下に既約ゴッパ符号を構成、復号のための必要最小限のコードが有ります。
コードは短いのですぐ読めると思います。


# 20230106

https://hal.inria.fr/inria-00073037/document

The Support Splitting Algorithm

とんでもない論文を見つけてしまった。さらに、

https://eprint.iacr.org/2022/514.pdf

この攻撃法を受けてパラメータを
$n=6688,q=8192,k=256$
に変更済み。

符号長が有限体のサイズに等しいときに攻撃できる（Goppa多項式が復元できる）とのことなので、短縮符号は安全らしい。

プログラムを全部非公開にするか迷ったけど、しばらく様子を見る。

Rustで暗号解読もやりたかったけど、今日選んだ自分の服の見立てが甘かったら、今設定している目標を変更しようと思っていたのだが。
まんざら間違っているとも言えなかったので、現在のテーマを進めることに。

グレブナー基底がわかったら２つもテーマに応用できるし、割り算ができるようになったことでさらに興味が出そう。

経済の矛盾。
商品価値？等価交換なんてそんなものは存在しない。
そんなことが事実だとしたら、利益はどこから出てくるのだ？

# 20230104

vc3000との結合テストが終わり、ほぼ完成しました。
vc3000はベクトル型のテストコードですが、今のところ無事に動いています。
多項式の表現をベクトル型に統一することで計算速度が飛躍的に早くなりました。
何をしているかというと秘密鍵である既約多項式を一発で出す計算です。

あとはわかりやすいようにコメントつけたり、変な関数名や変数名を直してバグを見つけて直したり、より簡潔なコードに仕上げて完成度を上げていきます。

# 20230103

とりあえずいろんな計算方法がごたまぜになったGoppaは引っ込めて、整理された形で始められるこのリポジトリを優先的に改造していくことにした。

もう一つは、二次元符号の使いみちだ。
これはまだはっきりしない。
年末は頭がぼーっとして考えようとすると眠くなるなど、なんだか調子が悪かった。

そんな中で新年も無事に迎えられて、とりあえずめでたい。

なんだか正月らしくも年末らしくもなく、ただ去年からつながった普通の一週間が過ぎただけのような感じのする今年の始まり。

健康不安があるので今年は御札を変えた。
健康が一番大事。
暗号を作れるのも健康あってこそ。

コロナから始まって、ウクライナの戦争、日本民主主義の崩壊の予感。

そんな恐怖が年末年始の私を取り巻いていた。
テレビでは新しい戦前が始まると言っている。

国家権力というものは基本的に邪悪なものだから、国民は常に警戒しているべきだというのが私の持論だが、みんな国のやることを警戒することもなく受け入れてしまっている。

こういう無防備さはとても危険だと思う。
そういう私も、暗号を作ったらダメだとか、作ったから逮捕だ、と言われたら抵抗のしようがない。
正当防衛を認めないのと同じことだ。

何が悪いことかを決めるのは自分じゃないからだ。

かつてウィニーの事件があったように、グレーゾーン技術ははっきりと線引きしてほしい。

ところで今日は、今年の年越しに、年越しそばを食べるのを忘れた事を思い出す一日であった。


# 20221231

とりあえず今一番完成形に近いのがこのリポジトリなので、るばーとさんと一緒に完成にまで持っていきたいです。

McElieceのRustバージョンがありますが、Rustのニーダーライターはまだないようなのでそれもやるかもしれないです。

### 国が民主主義に反することをしていても国民には対抗手段がないというのがとても不安です。

だって今日は平和でも明日は平和じゃないかもしれない。
それが一番不安です。

# 20220911

cruptanalysis of te niederreiter public key scheme based on grs subcodes

security bounds for the design of code-based cryptosystem

selecting parameters for secure mceliece-based cryptosystems

efficient implementation of the mceliece cryptosystem

attaking and defending the mceliece cryptosystem 

understanding binary-goppa decoding


# 20220826

これは暗号じゃないけど、暗号っていうのは隠れて使うものだという意見には疑問を感じる。
