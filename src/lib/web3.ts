import { http, createConfig } from 'wagmi'
import { mainnet, sepolia } from 'wagmi/chains'
import { injected } from 'wagmi/connectors'

// Story Network Testnet (Aenid)
export const storyTestnet = {
  id: 1315,
  name: 'Story Aenid',
  nativeCurrency: {
    decimals: 18,
    name: 'IP',
    symbol: 'IP',
  },
  rpcUrls: {
    default: { http: ['https://aeneid.storyrpc.io'] },
  },
  blockExplorers: {
    default: { name: 'Story Explorer', url: 'https://aeneid.storyscan.xyz' },
  },
  testnet: true,
} as const

export const config = createConfig({
  chains: [storyTestnet, mainnet, sepolia],
  connectors: [
    injected(),
  ],
  transports: {
    [storyTestnet.id]: http(),
    [mainnet.id]: http(),
    [sepolia.id]: http(),
  },
})
