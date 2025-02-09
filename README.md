# subspace
这是一个关于论文中描述的 MATLAB 算法的 GitHub 仓库。本仓库包含了论文中实验部分所使用的 MATLAB 代码。
## 目录
- [算法简介](#算法简介)
- [开发前的配置要求](#开发前的配置要求)
- [文件目录说明](#文件目录说明)
- [版本控制](#版本控制)
- [作者](#作者)
## 算法简介

高维线性逆问题（HDLIP）是不适定的，计算要求很高；然而，它们的解通常是稀疏的。基于具有非光滑重尾马蹄先验的贝叶斯框架，我们建立了一个子空间Gibbs采样方法来模拟两种高维线性逆问题(即变量选择和图像逆问题)的复杂后验。

## 开发前的配置要求

MATLAB R2018a及以上

## 安装

1. 克隆或下载本仓库。

2. 打开 MATLAB。

3. 导航到仓库的目录。

## 使用说明

```matlab
%Load a dataset:
[x1,x2,time] = bhs_compare(A, y_noisy,nsamples,burnin, thin);
```

### 文件目录说明

```
filetree 
├── ARCHITECTURE.md
├── LICENSE.txt
├── README.md
│  ├── /variable selection problems/
│  │  ├── bhs_compare.m
│  │  ├── test.m
│  │  └── dataset.mat
└── /util/

```


### 版本控制

该项目使用Git进行版本管理。您可以在repository参看当前可用版本。

### 作者

m202410737@xs.ustb.edu.cn
ijwmip@ustb.edu.cn
