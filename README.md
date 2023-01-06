# berlekamp-small

# 20230106

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
