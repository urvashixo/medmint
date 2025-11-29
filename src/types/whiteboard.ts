export interface Point {
  x: number
  y: number
  pressure?: number
}

export interface DrawingElement {
  id: string
  type: 'pen' | 'rectangle' | 'circle' | 'line' | 'text' | 'arrow' | 'diamond'
  points: Point[]
  color: string
  strokeWidth: number
  userId: string
  userName: string
  timestamp: number
  isDeleted?: boolean
  text?: string
  fontSize?: number
  roughness?: number
  fill?: string
  opacity?: number
  seed?: number
  strokeStyle?: 'solid' | 'dashed' | 'dotted'
  fillStyle?: 'solid' | 'hachure' | 'cross-hatch' | 'dots'
  angle?: number
  width?: number
  height?: number
  x?: number
  y?: number
}

export interface WhiteboardData {
  elements: DrawingElement[]
  version: number
  appState?: {
    viewBackgroundColor: string
    gridSize: number
    zoom: number
    scrollX: number
    scrollY: number
  }
}

export interface Cursor {
  userId: string
  userName: string
  x: number
  y: number
  color: string
  timestamp: number
}

export interface Tool {
  type: 'select' | 'pen' | 'rectangle' | 'circle' | 'line' | 'text' | 'eraser' | 'arrow' | 'diamond'
  color: string
  strokeWidth: number
  roughness: number
  fill: string
  opacity: number
  strokeStyle: 'solid' | 'dashed' | 'dotted'
  fillStyle: 'solid' | 'hachure' | 'cross-hatch' | 'dots'
}

export interface ViewState {
  zoom: number
  offsetX: number
  offsetY: number
}

export interface BoundingBox {
  x: number
  y: number
  width: number
  height: number
}