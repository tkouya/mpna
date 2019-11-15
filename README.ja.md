多倍長精度数値計算：サンプルプログラム
============================================================

2019-11-19 (Fri) 幸谷 智紀
---------------------------------

　以下のファイルは「[多倍長精度数値計算](https://www.morikita.co.jp/books/book/3393)」(森北出版)で解説しているプログラムです。詳細については本文を参照して下さい。

　本文及び演習問題で使われているプログラムはLinux環境下で実行を確認したものです。Windows 10環境下でも同等の環境を整えることができます。
 
☆コンパイル条件
-----------------------------

　本プログラムはLinuxソフトウェア開発環境下でコンパイル＆実行可能であることを下記の環境で確認しております。

・GCC 7.4.0 と  [GNU MP](https://gmplib.org/)(GMP), [MPFR](https://www.mpfr.org/), [QD](https://www.davidhbailey.com/dhbsoftware/), [MPFI](http://perso.ens-lyon.fr/nathalie.revol/software.html), [ARB](https://github.com/fredrik-johansson/arb), [mpreal.h](https://bitbucket.org/advanpix/mpreal/src/default/) がインストールされたUbuntu 18.04.01

・Intel compiler 13.1.3(GCC 4.4.7 compatible)と [GNU MP](https://gmplib.org/)(GMP), [MPFR](https://www.mpfr.org/), [QD](https://www.davidhbailey.com/dhbsoftware/), [MPFI](http://perso.ens-lyon.fr/nathalie.revol/software.html), [ARB](https://github.com/fredrik-johansson/arb), [mpreal.h](https://bitbucket.org/advanpix/mpreal/src/default/) がインストールされたCentOS 6.5


　Linux, Windows環境下での本プログラムのコンパイル＆実行は，上記のソフトウェア環境が整っているCUIで行って下さい。それ以外の環境下での諸問題については確認ができませんので，お答えすることも不可能です。


☆コンパイル方法
-----------------------------

1. Intel C/C++ Compilerの場合はicc.incを，GCCの場合はgcc.incを環境に合わせて修正し，それぞれのコンパイラが適切に動作するよう環境設定を行い，Makefileが読み込むファイルを設定
2. make でコンパイル  
3. make clean でobjectファイル，実行ファイルが消去される  
